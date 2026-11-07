% Independent MATLAB Antenna Toolbox reference for the coated-wire bend study.
%
% Keep this file as the stable, user-approved MATLAB entry point.  Runtime
% choices belong in tests/coated-wire/m7-bend-matlab-config.json so repeated
% invocations do not require broader permission to execute arbitrary MATLAB.

root = fileparts(fileparts(mfilename('fullpath')));
config_path = fullfile(root, 'tests', 'coated-wire', ...
    'm7-bend-matlab-config.json');
config = jsondecode(fileread(config_path));
if isfield(config, 'output_filename')
    output_filename = config.output_filename;
else
    output_filename = 'm7-bend-matlab-result.json';
end
output_path = fullfile(root, 'tests', 'coated-wire', output_filename);

result = struct;
result.schema = 'nec2c-coated-wire-m7-bend-matlab';
result.matlab_version = version;
result.operation = config.operation;
result.generated_utc = char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));

switch config.operation
    case 'inventory'
        result.antenna_toolbox = ver('antenna');
        result.custom_antenna_help = evalc('help customAntenna');
        result.custom_geometry_help = evalc('help customAntennaGeometry');
        result.cylinder_help = evalc('help shape.Cylinder');
        result.create_feed_help = evalc('help customAntenna/createFeed');
        result.mesh_help = evalc('help customAntenna/mesh');
        result.rotate_help = evalc('help shape.Cylinder/rotate');
    case 'smoke'
        if strcmp(config.geometry, 'straight')
            points = [0 0 -0.237; 0 0 0; 0 0 0.237];
        elseif strcmp(config.geometry, 'sharp-u')
            points = [-0.050 0 0.187; -0.050 0 0; 0 0 0; ...
                0.050 0 0; 0.050 0 0.187];
        else
            error('Unknown smoke-test geometry: %s', config.geometry);
        end
        body = polyline_shape(points, config.conductor_radius_m);
        ant = customAntenna(Shape=body);
        createFeed(ant, [0 0 0], 16);
        mesh_data = mesh(ant, MaxEdgeLength=config.maximum_edge_m);
        result.feed_location_m = ant.FeedLocation;
        result.geometry = config.geometry;
        result.mesh = struct('triangles', size(mesh_data.Triangles, 2), ...
            'points', size(mesh_data.Points, 2));
        frequencies = config.frequency_mhz(:) * 1e6;
        z = impedance(ant, frequencies);
        result.frequency_mhz = config.frequency_mhz(:);
        result.impedance_ohm = [real(z(:)), imag(z(:))];
    case 'study'
        frequencies = config.frequency_mhz(:) * 1e6;
        cases = cell(numel(config.cases), 1);
        for case_index = 1:numel(config.cases)
            source_case = config.cases(case_index);
            fprintf('MATLAB bend case %d/%d: %s\n', case_index, ...
                numel(config.cases), source_case.name);
            item = struct;
            item.name = source_case.name;
            item.points_m = source_case.points_m;
            started = tic;
            try
                body = polyline_shape(source_case.points_m, ...
                    config.conductor_radius_m);
                ant = customAntenna(Shape=body);
                createFeed(ant, [0 0 0], 16);
                mesh_data = mesh(ant, ...
                    MaxEdgeLength=config.maximum_edge_m);
                z = impedance(ant, frequencies);
                item.status = 'ok';
                item.mesh = struct( ...
                    'triangles', size(mesh_data.Triangles, 2), ...
                    'points', size(mesh_data.Points, 2));
                item.frequency_mhz = config.frequency_mhz(:);
                item.impedance_ohm = [real(z(:)), imag(z(:))];
                [f0, r0] = find_resonance(config.frequency_mhz(:), z(:));
                item.resonance = struct('frequency_mhz', f0, ...
                    'resistance_ohm', r0);
            catch exception
                item.status = 'error';
                item.error_identifier = exception.identifier;
                item.error_message = exception.message;
            end
            item.wall_seconds = toc(started);
            cases{case_index} = item;
        end
        result.cases = [cases{:}];
    otherwise
        error('Unknown coated-wire bend MATLAB operation: %s', ...
            config.operation);
end

text = jsonencode(result, PrettyPrint=true);
fid = fopen(output_path, 'w');
assert(fid >= 0, 'Could not open MATLAB result path');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', text);
fprintf('Wrote %s\n', output_path);

function body = polyline_shape(points, radius)
    % A watertight, planar, mitered tube avoids CAD Boolean ambiguity where
    % perpendicular closed cylinders share an end face.
    sides = 16;
    count = size(points, 1);
    vertices = zeros(count*sides+2, 3);
    y_normal = [0 1 0];
    for k = 1:count
        if k == 1
            tangent = points(2, :)-points(1, :);
        elseif k == count
            tangent = points(count, :)-points(count-1, :);
        else
            incoming = points(k, :)-points(k-1, :);
            outgoing = points(k+1, :)-points(k, :);
            tangent = incoming/norm(incoming)+outgoing/norm(outgoing);
        end
        tangent = tangent/norm(tangent);
        radial = cross(tangent, y_normal);
        radial = radial/norm(radial);
        for j = 1:sides
            angle = 2*pi*(j-1)/sides;
            vertices((k-1)*sides+j, :) = points(k, :) + radius*( ...
                cos(angle)*y_normal + sin(angle)*radial);
        end
    end
    first_center = count*sides+1;
    last_center = count*sides+2;
    vertices(first_center, :) = points(1, :);
    vertices(last_center, :) = points(end, :);
    faces = zeros(2*(count-1)*sides+2*sides, 3);
    f = 0;
    for k = 1:count-1
        for j = 1:sides
            next = mod(j, sides)+1;
            a = (k-1)*sides+j;
            b = (k-1)*sides+next;
            c = k*sides+j;
            d = k*sides+next;
            f = f+1; faces(f, :) = [a c b];
            f = f+1; faces(f, :) = [b c d];
        end
    end
    for j = 1:sides
        next = mod(j, sides)+1;
        f = f+1; faces(f, :) = [first_center next j];
        a = (count-1)*sides+j;
        b = (count-1)*sides+next;
        f = f+1; faces(f, :) = [last_center a b];
    end
    body = shape.Custom3D(triangulation(faces, vertices));
end

function [frequency, resistance] = find_resonance(frequencies, impedance)
    frequency = NaN;
    resistance = NaN;
    for k = 1:numel(frequencies)-1
        if imag(impedance(k)) == 0 || ...
                imag(impedance(k))*imag(impedance(k+1)) < 0
            if imag(impedance(k)) == 0
                part = 0;
            else
                part = -imag(impedance(k))/(imag(impedance(k+1))- ...
                    imag(impedance(k)));
            end
            frequency = frequencies(k) + part*(frequencies(k+1)- ...
                frequencies(k));
            resistance = real(impedance(k)) + part*( ...
                real(impedance(k+1))-real(impedance(k)));
            return
        end
    end
end
