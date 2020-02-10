# Directories where we can build something.
SUBDIRS          = doc lib-src src

# Install the examples into /usr/share/doc/cddlib/examples*
nobase_doc_DATA = $(srcdir)/examples/* $(srcdir)/examples-ine/* $(srcdir)/examples-ine3d/* $(srcdir)/examples-ext/*

# Additional bits that should be part of the released tarball.
EXTRA_DIST = $(nobase_doc_DATA) $(srcdir)/README.core2processor $(srcdir)/bootstrap

ACLOCAL_AMFLAGS = -I m4

# All localdebug should be disabled before we package a releasable tarball
distcheck-hook:
	! find $(srcdir) -type f -print0 | xargs -0 $(EGREP) 'localdebug *= *dd_TRUE'

pkgconfigdir       = $(libdir)/pkgconfig
pkgconfig_DATA     = cddlib.pc
