# Directories where we need to build something.
if GMP
    SUBDIRS      = lib-src src lib-src-gmp src-gmp
else
    # Do not build GMP libraries if GMP is not available.
    SUBDIRS      = lib-src src
endif

# Additional directories that should be part of the released tarball.
EXTRA_DIST   = $(srcdir)/examples/*.* $(srcdir)/examples-ine/*.* \
               $(srcdir)/examples-ine3d/*.* $(srcdir)/examples-ext/*.*

ACLOCAL_AMFLAGS = -I m4
