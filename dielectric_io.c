#include "nec2c.h"
#include "dielectric.h"

#include <ctype.h>
#include <limits.h>
#include <stdint.h>

#define TOKEN_LIMIT 32

static void set_diagnostic(dielectric_diagnostic *diagnostic,
                           const char *resource, size_t line,
                           const char *record, const char *value,
                           const char *message)
{
  snprintf(diagnostic->resource, sizeof(diagnostic->resource), "%s",
           resource != NULL ? resource : "<DF card>");
  snprintf(diagnostic->record, sizeof(diagnostic->record), "%s",
           record != NULL ? record : "UNKNOWN");
  snprintf(diagnostic->value, sizeof(diagnostic->value), "%s",
           value != NULL ? value : "");
  snprintf(diagnostic->message, sizeof(diagnostic->message), "%s", message);
  diagnostic->line = line;
}

static int keyword_equal(const char *left, const char *right)
{
  while((*left != '\0') && (*right != '\0'))
  {
    if(toupper((unsigned char)*left) != toupper((unsigned char)*right))
      return 0;
    left++;
    right++;
  }
  return (*left == '\0') && (*right == '\0');
}

static char *copy_string(const char *text)
{
  size_t length = strlen(text) + 1;
  char *copy = malloc(length);
  if(copy != NULL)
    memcpy(copy, text, length);
  return copy;
}

static int parse_positive_int(const char *text, int *value)
{
  char *end;
  long parsed;
  errno = 0;
  parsed = strtol(text, &end, 10);
  if((errno != 0) || (*text == '\0') || (*end != '\0') ||
     (parsed <= 0) || (parsed > INT_MAX))
    return 0;
  *value = (int)parsed;
  return 1;
}

static int parse_positive_size(const char *text, size_t *value)
{
  char *end;
  unsigned long long parsed;
  errno = 0;
  if((*text == '\0') || (*text == '-') || (*text == '+'))
    return 0;
  parsed = strtoull(text, &end, 10);
  if((errno != 0) || (*end != '\0') || (parsed == 0) ||
     (parsed > (unsigned long long)SIZE_MAX))
    return 0;
  *value = (size_t)parsed;
  return 1;
}

static int parse_finite_double(const char *text, double *value)
{
  char *end;
  errno = 0;
  *value = strtod(text, &end);
  return (errno == 0) && (*text != '\0') && (*end == '\0') && isfinite(*value);
}

static int tokenize(char *line, char **tokens, int *unterminated_quote)
{
  int count = 0;
  char *cursor = line;
  *unterminated_quote = 0;

  while(*cursor != '\0')
  {
    while(isspace((unsigned char)*cursor))
      cursor++;
    if((*cursor == '\0') || (*cursor == '#'))
      break;
    if(count == TOKEN_LIMIT)
      return -1;
    if(*cursor == '"')
    {
      tokens[count++] = ++cursor;
      while((*cursor != '\0') && (*cursor != '"'))
        cursor++;
      if(*cursor != '"')
      {
        *unterminated_quote = 1;
        return count;
      }
      *cursor++ = '\0';
      if((*cursor != '\0') && (*cursor != '#') &&
         !isspace((unsigned char)*cursor))
        return -1;
    }
    else
    {
      tokens[count++] = cursor;
      while((*cursor != '\0') && (*cursor != '#') &&
            !isspace((unsigned char)*cursor))
        cursor++;
      if(*cursor == '#')
      {
        *cursor = '\0';
        break;
      }
      if(*cursor != '\0')
        *cursor++ = '\0';
    }
  }
  return count;
}

static int parse_df_path(const char *field, char *path, size_t path_size,
                         dielectric_diagnostic *diagnostic)
{
  const char *cursor = field;
  const char *start;
  size_t length;
  while(isspace((unsigned char)*cursor))
    cursor++;
  if(*cursor == '"')
  {
    cursor++;
    start = cursor;
    while((*cursor != '\0') && (*cursor != '"'))
      cursor++;
    if(*cursor != '"')
    {
      set_diagnostic(diagnostic, NULL, 0, "DF", field,
                     "unterminated quoted resource path");
      return 0;
    }
    length = (size_t)(cursor - start);
    cursor++;
  }
  else
  {
    start = cursor;
    while((*cursor != '\0') && !isspace((unsigned char)*cursor))
      cursor++;
    length = (size_t)(cursor - start);
  }
  while(isspace((unsigned char)*cursor))
    cursor++;
  if((length == 0) || (*cursor != '\0'))
  {
    set_diagnostic(diagnostic, NULL, 0, "DF", field,
                   length == 0 ? "resource path is missing" :
                                 "unexpected content after resource path");
    return 0;
  }
  if(length >= path_size)
  {
    set_diagnostic(diagnostic, NULL, 0, "DF", field,
                   "resource path exceeds the supported path limit");
    return 0;
  }
  memcpy(path, start, length);
  path[length] = '\0';
  return 1;
}

static int path_is_absolute(const char *path)
{
  if((path[0] == '/') || (path[0] == '\\'))
    return 1;
  return isalpha((unsigned char)path[0]) && (path[1] == ':') &&
         ((path[2] == '/') || (path[2] == '\\'));
}

static int resolve_path(const char *deck_path, const char *resource_path,
                        char *resolved, size_t resolved_size)
{
  const char *slash;
  size_t directory_length;
  if(path_is_absolute(resource_path))
    return snprintf(resolved, resolved_size, "%s", resource_path) <
           (int)resolved_size;
  slash = strrchr(deck_path, '/');
  if(slash == NULL)
    slash = strrchr(deck_path, '\\');
  directory_length = slash == NULL ? 0 : (size_t)(slash - deck_path + 1);
  if(directory_length + strlen(resource_path) + 1 > resolved_size)
    return 0;
  memcpy(resolved, deck_path, directory_length);
  strcpy(resolved + directory_length, resource_path);
  return 1;
}

static int append_material(dielectric_model *model, dielectric_material value)
{
  dielectric_material *items = realloc(model->materials,
      (model->material_count + 1) * sizeof(*items));
  if(items == NULL)
    return 0;
  model->materials = items;
  model->materials[model->material_count++] = value;
  return 1;
}

static int append_region(dielectric_model *model, dielectric_region value)
{
  dielectric_region *items = realloc(model->regions,
      (model->region_count + 1) * sizeof(*items));
  if(items == NULL)
    return 0;
  model->regions = items;
  model->regions[model->region_count++] = value;
  return 1;
}

static int append_box(dielectric_model *model, dielectric_box value)
{
  dielectric_box *items = realloc(model->boxes,
      (model->box_count + 1) * sizeof(*items));
  if(items == NULL)
    return 0;
  model->boxes = items;
  model->boxes[model->box_count++] = value;
  return 1;
}

static int material_exists(const dielectric_model *model, int id)
{
  size_t index;
  for(index = 0; index < model->material_count; index++)
    if(model->materials[index].id == id)
      return 1;
  return 0;
}

static int region_exists(const dielectric_model *model, int id)
{
  size_t index;
  for(index = 0; index < model->region_count; index++)
    if(model->regions[index].id == id)
      return 1;
  return 0;
}

static int compare_materials(const void *left, const void *right)
{
  const dielectric_material *a = left;
  const dielectric_material *b = right;
  return (a->id > b->id) - (a->id < b->id);
}

static int compare_regions(const void *left, const void *right)
{
  const dielectric_region *a = left;
  const dielectric_region *b = right;
  return (a->id > b->id) - (a->id < b->id);
}

static int compare_boxes(const void *left, const void *right)
{
  const dielectric_box *a = left;
  const dielectric_box *b = right;
  int axis;
  if(a->region_id != b->region_id)
    return (a->region_id > b->region_id) - (a->region_id < b->region_id);
  for(axis = 0; axis < 3; axis++)
  {
    if(a->minimum[axis] != b->minimum[axis])
      return a->minimum[axis] > b->minimum[axis] ? 1 : -1;
    if(a->maximum[axis] != b->maximum[axis])
      return a->maximum[axis] > b->maximum[axis] ? 1 : -1;
  }
  return 0;
}

void dielectric_model_init(dielectric_model *model)
{
  memset(model, 0, sizeof(*model));
  model->max_kh = DIELECTRIC_DEFAULT_MAX_KH;
  model->max_cells = DIELECTRIC_DEFAULT_MAX_CELLS;
  model->max_aspect = DIELECTRIC_DEFAULT_MAX_ASPECT;
}

void dielectric_model_free(dielectric_model *model)
{
  free(model->resource_path);
  free(model->materials);
  free(model->regions);
  free(model->boxes);
  dielectric_model_init(model);
}

static int parse_resource(dielectric_model *model, FILE *resource,
                          dielectric_diagnostic *diagnostic)
{
  char line[DIELECTRIC_RESOURCE_LINE_MAX + 2];
  size_t line_number = 0;
  int header_seen = 0;
  int end_seen = 0;

  while(fgets(line, sizeof(line), resource) != NULL)
  {
    char *tokens[TOKEN_LIMIT];
    int token_count;
    int unterminated_quote;
    size_t length;
    line_number++;
    length = strlen(line);
    if((length > 0) && (line[length - 1] != '\n') && !feof(resource))
    {
      int character;
      while(((character = fgetc(resource)) != '\n') && (character != EOF));
      set_diagnostic(diagnostic, model->resource_path, line_number, "LINE", "",
                     "line exceeds 1024 characters");
      return 0;
    }
    token_count = tokenize(line, tokens, &unterminated_quote);
    if(unterminated_quote)
    {
      set_diagnostic(diagnostic, model->resource_path, line_number, "LINE", "\"",
                     "unterminated quoted token");
      return 0;
    }
    if(token_count < 0)
    {
      set_diagnostic(diagnostic, model->resource_path, line_number, "LINE", "",
                     "too many tokens");
      return 0;
    }
    if(token_count == 0)
      continue;
    if(end_seen)
    {
      set_diagnostic(diagnostic, model->resource_path, line_number, tokens[0],
                     tokens[0], "content is not permitted after END");
      return 0;
    }
    if(!header_seen)
    {
      int version;
      if(!keyword_equal(tokens[0], "NEC2C-DIELECTRIC") || (token_count != 2))
      {
        set_diagnostic(diagnostic, model->resource_path, line_number, tokens[0],
                       tokens[0], "first record must be NEC2C-DIELECTRIC <version>");
        return 0;
      }
      if(!parse_positive_int(tokens[1], &version) || (version != 1))
      {
        set_diagnostic(diagnostic, model->resource_path, line_number,
                       "NEC2C-DIELECTRIC", tokens[1],
                       "unsupported format version; expected 1");
        return 0;
      }
      model->format_version = version;
      header_seen = 1;
      continue;
    }
    if(keyword_equal(tokens[0], "MATERIAL"))
    {
      dielectric_material material;
      if((token_count != 4) || !keyword_equal(tokens[2], "ER"))
      {
        set_diagnostic(diagnostic, model->resource_path, line_number, "MATERIAL",
                       tokens[0], "expected MATERIAL <id> ER <epsilon-r>");
        return 0;
      }
      if(!parse_positive_int(tokens[1], &material.id))
      {
        set_diagnostic(diagnostic, model->resource_path, line_number, "MATERIAL",
                       tokens[1], "material ID must be a positive decimal integer");
        return 0;
      }
      if(material_exists(model, material.id))
      {
        set_diagnostic(diagnostic, model->resource_path, line_number, "MATERIAL",
                       tokens[1], "duplicate material ID");
        return 0;
      }
      if(!parse_finite_double(tokens[3], &material.epsilon_r) ||
         (material.epsilon_r < 1.0))
      {
        set_diagnostic(diagnostic, model->resource_path, line_number, "MATERIAL",
                       tokens[3], "ER must be finite and greater than or equal to 1");
        return 0;
      }
      material.source_line = line_number;
      if(!append_material(model, material))
        goto allocation_failure;
    }
    else if(keyword_equal(tokens[0], "REGION"))
    {
      dielectric_region region;
      if((token_count != 4) || !keyword_equal(tokens[2], "MATERIAL"))
      {
        set_diagnostic(diagnostic, model->resource_path, line_number, "REGION",
                       tokens[0], "expected REGION <id> MATERIAL <material-id>");
        return 0;
      }
      if(!parse_positive_int(tokens[1], &region.id) ||
         !parse_positive_int(tokens[3], &region.material_id))
      {
        set_diagnostic(diagnostic, model->resource_path, line_number, "REGION",
                       !parse_positive_int(tokens[1], &region.id) ? tokens[1] : tokens[3],
                       "region and material IDs must be positive decimal integers");
        return 0;
      }
      if(region_exists(model, region.id))
      {
        set_diagnostic(diagnostic, model->resource_path, line_number, "REGION",
                       tokens[1], "duplicate region ID");
        return 0;
      }
      region.source_line = line_number;
      if(!append_region(model, region))
        goto allocation_failure;
    }
    else if(keyword_equal(tokens[0], "BOX"))
    {
      dielectric_box box;
      int axis;
      if(token_count != 8)
      {
        set_diagnostic(diagnostic, model->resource_path, line_number, "BOX",
                       tokens[0], "expected BOX <region-id> xmin ymin zmin xmax ymax zmax");
        return 0;
      }
      if(!parse_positive_int(tokens[1], &box.region_id))
      {
        set_diagnostic(diagnostic, model->resource_path, line_number, "BOX",
                       tokens[1], "region ID must be a positive decimal integer");
        return 0;
      }
      for(axis = 0; axis < 3; axis++)
      {
        if(!parse_finite_double(tokens[2 + axis], &box.minimum[axis]))
        {
          set_diagnostic(diagnostic, model->resource_path, line_number, "BOX",
                         tokens[2 + axis], "box bounds must be finite numbers");
          return 0;
        }
        if(!parse_finite_double(tokens[5 + axis], &box.maximum[axis]))
        {
          set_diagnostic(diagnostic, model->resource_path, line_number, "BOX",
                         tokens[5 + axis], "box bounds must be finite numbers");
          return 0;
        }
        if(box.minimum[axis] >= box.maximum[axis])
        {
          set_diagnostic(diagnostic, model->resource_path, line_number, "BOX",
                         tokens[2 + axis], "each minimum bound must be less than its maximum");
          return 0;
        }
      }
      box.source_line = line_number;
      if(!append_box(model, box))
        goto allocation_failure;
    }
    else if(keyword_equal(tokens[0], "MESH"))
    {
      int index;
      int seen_kh = 0, seen_cells = 0, seen_aspect = 0;
      if(model->mesh_record_present)
      {
        set_diagnostic(diagnostic, model->resource_path, line_number, "MESH",
                       tokens[0], "duplicate MESH record");
        return 0;
      }
      if(((token_count - 1) % 2) != 0)
      {
        set_diagnostic(diagnostic, model->resource_path, line_number, "MESH",
                       tokens[token_count - 1], "mesh options require keyword/value pairs");
        return 0;
      }
      for(index = 1; index < token_count; index += 2)
      {
        if(keyword_equal(tokens[index], "MAX_KH"))
        {
          if(seen_kh++)
            goto duplicate_mesh_option;
          if(!parse_finite_double(tokens[index + 1], &model->max_kh) ||
             (model->max_kh <= 0.0))
            goto invalid_mesh_value;
        }
        else if(keyword_equal(tokens[index], "MAX_CELLS"))
        {
          if(seen_cells++)
            goto duplicate_mesh_option;
          if(!parse_positive_size(tokens[index + 1], &model->max_cells))
            goto invalid_mesh_value;
        }
        else if(keyword_equal(tokens[index], "MAX_ASPECT"))
        {
          if(seen_aspect++)
            goto duplicate_mesh_option;
          if(!parse_finite_double(tokens[index + 1], &model->max_aspect) ||
             (model->max_aspect < 1.0))
            goto invalid_mesh_value;
        }
        else
        {
          set_diagnostic(diagnostic, model->resource_path, line_number, "MESH",
                         tokens[index], "unknown mesh option");
          return 0;
        }
        continue;
duplicate_mesh_option:
        set_diagnostic(diagnostic, model->resource_path, line_number, "MESH",
                       tokens[index], "duplicate mesh option");
        return 0;
invalid_mesh_value:
        set_diagnostic(diagnostic, model->resource_path, line_number, "MESH",
                       tokens[index + 1], "invalid mesh option value");
        return 0;
      }
      model->mesh_record_present = 1;
    }
    else if(keyword_equal(tokens[0], "END"))
    {
      if(token_count != 1)
      {
        set_diagnostic(diagnostic, model->resource_path, line_number, "END",
                       tokens[1], "END does not accept arguments");
        return 0;
      }
      end_seen = 1;
    }
    else
    {
      set_diagnostic(diagnostic, model->resource_path, line_number, tokens[0],
                     tokens[0], "unknown record keyword");
      return 0;
    }
  }

  if(ferror(resource))
  {
    set_diagnostic(diagnostic, model->resource_path, line_number, "FILE", "",
                   "I/O error while reading resource");
    return 0;
  }
  if(!header_seen || !end_seen || (model->material_count == 0) ||
     (model->region_count == 0) || (model->box_count == 0))
  {
    const char *missing = !header_seen ? "NEC2C-DIELECTRIC" :
                          !end_seen ? "END" :
                          model->material_count == 0 ? "MATERIAL" :
                          model->region_count == 0 ? "REGION" : "BOX";
    set_diagnostic(diagnostic, model->resource_path, line_number, "FILE", missing,
                   "required record is missing");
    return 0;
  }
  {
    size_t index;
    for(index = 0; index < model->region_count; index++)
      if(!material_exists(model, model->regions[index].material_id))
      {
        char value[32];
        snprintf(value, sizeof(value), "%d", model->regions[index].material_id);
        set_diagnostic(diagnostic, model->resource_path,
                       model->regions[index].source_line, "REGION", value,
                       "referenced material ID is not declared");
        return 0;
      }
    for(index = 0; index < model->box_count; index++)
      if(!region_exists(model, model->boxes[index].region_id))
      {
        char value[32];
        snprintf(value, sizeof(value), "%d", model->boxes[index].region_id);
        set_diagnostic(diagnostic, model->resource_path,
                       model->boxes[index].source_line, "BOX", value,
                       "referenced region ID is not declared");
        return 0;
      }
  }
  qsort(model->materials, model->material_count, sizeof(*model->materials), compare_materials);
  qsort(model->regions, model->region_count, sizeof(*model->regions), compare_regions);
  qsort(model->boxes, model->box_count, sizeof(*model->boxes), compare_boxes);
  return 1;

allocation_failure:
  set_diagnostic(diagnostic, model->resource_path, line_number, "FILE", "",
                 "memory allocation failed");
  return 0;
}

int dielectric_load(dielectric_model *model, const char *deck_path,
                    const char *df_field, dielectric_diagnostic *diagnostic)
{
  char resource_path[NEC2C_PATH_MAX * 4 + 1];
  char resolved[NEC2C_PATH_MAX * 4 + 1];
  FILE *resource;
  dielectric_model candidate;
  dielectric_model_init(&candidate);
  if(!parse_df_path(df_field, resource_path, sizeof(resource_path), diagnostic))
    return 0;
  if(!resolve_path(deck_path, resource_path, resolved, sizeof(resolved)))
  {
    set_diagnostic(diagnostic, NULL, 0, "DF", resource_path,
                   "resolved resource path exceeds the supported path limit");
    return 0;
  }
  candidate.resource_path = copy_string(resolved);
  candidate.path_is_absolute = path_is_absolute(resource_path);
  if(candidate.resource_path == NULL)
  {
    set_diagnostic(diagnostic, resolved, 0, "DF", resource_path,
                   "memory allocation failed");
    return 0;
  }
  resource = fopen(resolved, "r");
  if(resource == NULL)
  {
    set_diagnostic(diagnostic, resolved, 0, "DF", resource_path,
                   "cannot open dielectric resource");
    dielectric_model_free(&candidate);
    return 0;
  }
  if(!parse_resource(&candidate, resource, diagnostic))
  {
    fclose(resource);
    dielectric_model_free(&candidate);
    return 0;
  }
  fclose(resource);
  candidate.loaded = 1;
  dielectric_model_free(model);
  *model = candidate;
  return 1;
}

void dielectric_print_diagnostic(FILE *stream,
                                 const dielectric_diagnostic *diagnostic)
{
  fprintf(stream, "\n\n  DIELECTRIC RESOURCE ERROR: %s:%lu: %s: %s",
          diagnostic->resource, (unsigned long)diagnostic->line,
          diagnostic->record, diagnostic->message);
  if(diagnostic->value[0] != '\0')
    fprintf(stream, " [value: %s]", diagnostic->value);
}

void dielectric_print_summary(FILE *stream, const dielectric_model *model)
{
  size_t index;
  fprintf(stream, "\n\n\n"
          "                       -------- DIELECTRIC RESOURCE --------\n"
          "  FORMAT VERSION: %d\n"
          "  RESOLVED RESOURCE: %s\n"
          "  PATH PORTABILITY: %s\n"
          "  MATERIAL COUNT: %lu\n"
          "  REGION COUNT: %lu\n"
          "  PHYSICAL BOX COUNT: %lu\n"
          "\n  MATERIALS (ID, ER)",
          model->format_version, model->resource_path,
          model->path_is_absolute ? "ABSOLUTE (NONPORTABLE)" : "DECK-RELATIVE",
          (unsigned long)model->material_count,
          (unsigned long)model->region_count,
          (unsigned long)model->box_count);
  for(index = 0; index < model->material_count; index++)
    fprintf(stream, "\n  MATERIAL %d ER %.17g", model->materials[index].id,
            model->materials[index].epsilon_r);
  fprintf(stream, "\n\n  REGIONS (ID, MATERIAL)");
  for(index = 0; index < model->region_count; index++)
    fprintf(stream, "\n  REGION %d MATERIAL %d", model->regions[index].id,
            model->regions[index].material_id);
  fprintf(stream, "\n\n  PHYSICAL BOXES (REGION, XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX)");
  for(index = 0; index < model->box_count; index++)
    fprintf(stream, "\n  BOX %d %.17g %.17g %.17g %.17g %.17g %.17g",
            model->boxes[index].region_id,
            model->boxes[index].minimum[0], model->boxes[index].minimum[1],
            model->boxes[index].minimum[2], model->boxes[index].maximum[0],
            model->boxes[index].maximum[1], model->boxes[index].maximum[2]);
  fprintf(stream, "\n\n  MESH MAX_KH %.17g MAX_CELLS %lu MAX_ASPECT %.17g%s\n"
          "  STATUS: PARSED; READY FOR DIELECTRIC SOLVER",
          model->max_kh, (unsigned long)model->max_cells, model->max_aspect,
          model->mesh_record_present ? "" : " (DEFAULTS)");
}
