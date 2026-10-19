## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = subdir-objects

EXTRA_DIST = \
  autogen.sh \
  NEC2-bug.txt \
  docs/DIELECTRIC-RESOURCE-FORMAT.md \
  docs/VIE-CUBOID-KERNEL.md \
  docs/DIELECTRIC-SOLVER-M5.md \
  docs/M5-REFINEMENT-REPORT.md \
	  docs/DIELECTRIC-WIRE-M6.md \
	  docs/M6-OPENEMS-REPORT.md \
	  docs/DIELECTRIC-PATCH-M7.md \
  Input/EX1.nec \
  Input/EX4.nec \
  Input/EX5.nec \
  Input/EX7.nec \
  tests/README.md \
  tests/baseline.json \
  tests/dielectric/expected-summary.txt \
  tests/dielectric/expected-mesh.txt \
  tests/dielectric/resources/valid-model.diel \
  tests/dielectric/valid.nec \
  tests/dielectric/m5-cube.diel \
  tests/dielectric/m5-cube.nec \
	  tests/dielectric/m6-legacy.nec \
	  tests/dielectric/m6-vacuum.nec \
	  tests/dielectric/m6-vacuum.diel \
	  tests/dielectric/m6-cube.nec \
	  tests/dielectric/m6-cube.diel \
	  tests/dielectric/m6-cube-fine.nec \
	  tests/dielectric/m6-cube-fine.diel \
	  tests/dielectric/m6-cube-far.nec \
	  tests/dielectric/m6-cube-far.diel \
	  tests/dielectric/m6-vacuum-21.nec \
	  tests/dielectric/m6-cube-21.nec \
	  tests/dielectric/m6-n2.diel \
	  tests/dielectric/m6-n2.nec \
	  tests/dielectric/m6-n2-300.nec \
	  tests/dielectric/m6-n2-legacy.nec \
	  tests/dielectric/m6-n2-legacy-sweep.nec \
	  tests/dielectric/m7-patch-legacy.nec \
	  tests/dielectric/m7-patch-vacuum.nec \
	  tests/dielectric/m7-vacuum.diel \
	  tests/dielectric/m7-patch-cube.nec \
	  tests/dielectric/m7-patch-cube-coarse.nec \
	  tests/dielectric/m7-patch-cube-far.nec \
	  tests/dielectric/m7-patch-cube-er2.nec \
	  tests/dielectric/m7-cube.diel \
	  tests/dielectric/m7-cube-coarse.diel \
	  tests/dielectric/m7-cube-far.diel \
	  tests/dielectric/m7-cube-er2.diel \
	  tests/dielectric/m7-patch4-vacuum.nec \
	  tests/dielectric/m7-patch4-cube.nec \
	  tests/dielectric/m7-mixed-legacy.nec \
	  tests/dielectric/m7-mixed-vacuum.nec \
	  tests/dielectric/m7-mixed-vacuum.diel \
	  tests/dielectric/m7-mixed-cube.nec \
	  tests/dielectric/m7-mixed-cube.diel \
	  tests/dielectric/m7-mixed-cube-fine.nec \
	  tests/dielectric/m7-mixed-cube-fine.diel \
  tests/dielectric_parser.py \
  tests/dielectric_mesh.py \
  tests/dielectric_solver.py \
	  tests/dielectric_wire.py \
	  tests/dielectric_patch.py \
  tests/generate_vie_reference.py \
  tests/voxscatter-reference.csv \
	  tests/openems-m6-n2-15mm.csv \
	  tests/openems-m6-n2-10mm.csv \
	  tests/nec2c-m6-n2.csv \
	  tests/nec2c-m7-mixed.csv \
	  tests/openems-m7-mixed-15mm.csv \
	  tests/openems-m7-mixed-10mm.csv \
  tests/vie-reference.csv \
  tests/vie_kernel_test.c \
  tools/directfn_surface_oracle.cpp \
  tools/m5_refinement.py \
	  tools/openems_m6_n2.py \
	  tools/openems_m7_mixed.py \
  tools/voxscatter_m5_reference.m \
  tests/decks/nt-network.nec \
  tests/decks/patch-only.nec \
  tests/regression.py

TESTS = tests/regression.py tests/dielectric_parser.py tests/dielectric_mesh.py \
	tests/dielectric_solver.py \
	tests/dielectric_wire.py \
	tests/dielectric_patch.py \
	tests/vie_kernel_test
TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON)
TESTS_ENVIRONMENT = \
	NEC2C='$(abs_top_builddir)/nec2c'; export NEC2C; \
	TEST_SRCDIR='$(abs_top_srcdir)'; export TEST_SRCDIR;

.PHONY: regression
regression: nec2c
	$(TESTS_ENVIRONMENT) $(PYTHON) $(srcdir)/tests/regression.py

AM_CPPFLAGS = \
	-DPACKAGE_DATA_DIR=\""$(datadir)"\" \
	-DPACKAGE_LOCALE_DIR=\""$(prefix)/$(DATADIRNAME)/locale"\" \
	@PACKAGE_CFLAGS@

bin_PROGRAMS = nec2c

check_PROGRAMS = tests/vie_kernel_test
tests_vie_kernel_test_SOURCES = tests/vie_kernel_test.c vie_kernel.c vie_kernel.h
tests_vie_kernel_test_LDADD = -lm

nec2c_SOURCES = \
  calculations.c \
  dielectric_io.c \
  dielectric_mesh.c \
  dielectric_solver.c \
  dielectric_wire.c \
  dielectric.h \
  vie_kernel.c \
  vie_kernel.h \
  geometry.c \
  input.c \
  matrix.c \
  network.c \
  shared.c \
  fields.c \
  ground.c \
  main.c \
  misc.c \
  radiation.c \
  somnec.c \
  nec2c.h \
  shared.h

install-data-local:
	@$(NORMAL_INSTALL)
	if test -d $(srcdir)/pixmaps; then \
	  $(mkinstalldirs) $(DESTDIR)$(pkgdatadir)/pixmaps; \
	  for pixmap in $(srcdir)/pixmaps/*; do \
	    if test -f $$pixmap; then \
	      $(INSTALL_DATA) $$pixmap $(DESTDIR)$(pkgdatadir)/pixmaps; \
	    fi \
	  done \
	fi

dist-hook:
	if test -d pixmaps; then \
	  mkdir $(distdir)/pixmaps; \
	  for pixmap in pixmaps/*; do \
	    if test -f $$pixmap; then \
	      cp -p $$pixmap $(distdir)/pixmaps; \
	    fi \
	  done \
	fi

