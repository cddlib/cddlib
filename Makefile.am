# Do not build GMP libraries if GMP is not available.
if GMP
    GMP_SUBDIRS  = lib-src-gmp src-gmp
endif

# Directories where we can build something.
SUBDIRS          = doc lib-src src $(GMP_SUBDIRS)

# Additional directories that should be part of the released tarball.
EXTRA_DIST   = $(srcdir)/examples $(srcdir)/examples-ine $(srcdir)/examples-ine3d $(srcdir)/examples-ext

ACLOCAL_AMFLAGS = -I m4
