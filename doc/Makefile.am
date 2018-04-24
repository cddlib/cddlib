# The Makefile for doc/ is non-standard as we distribute a manual that is
# built from a .tex file (and not a .texi file.)
# We do not want to make latex & pdflatx & dvips requirements for building
# cddlib, so we ship the built manual as part of the distribution.

if LATEX
LATEX_SRC = $(srcdir)/cddlibman.tex $(srcdir)/cddlibman.bbl
# If these auxiliary files are up-to-date, then we'll most likely get the labels right.
LATEX_AUX = cddlibman.toc cddlibman.aux

# Building twice seems to be enough to get references right.
$(LATEX_AUX): ${LATEX_SRC}
	$(latex) $<
	$(latex) $<

cddlibman.dvi: ${LATEX_SRC} | ${LATEX_AUX}
	$(latex) $<

if PDF
cddlibman.pdf: ${LATEX_SRC} | ${LATEX_AUX}
	$(pdflatex) $<
endif

if PS
cddlibman.ps: cddlibman.dvi
	$(dvips) cddlibman.dvi
endif

if HTML
# Wipe the HTML output directory as latex2html's perl calls otherwise try to
# find images.pl in there without putting cddlibman into @INC.
cddlibman/cddlibman.html: ${LATEX_SRC} html.sty | ${LATEX_AUX}
	rm -rf cddlibman
	$(latex2html) $<
endif

endif

# Include the PDF, PS, and DVI documentation in the distributed tarball.  Note
# that this means that pdflatex & latex & dvips must be available on the system
# where make dist(check) is run.
EXTRA_DIST = $(srcdir)/cddlibman.pdf $(srcdir)/cddlibman.dvi $(srcdir)/cddlibman.ps \
             $(LATEX_SRC) $(srcdir)/html.sty

# Clean up files that are not in the tarball
CLEANFILES = cddlibman.toc cddlibman.aux cddlibman.log
# Clean up files that "make dist" generates
MAINTAINERCLEANFILES = cddlibman.pdf cddlibman.ps cddlibman.dvi
# Change this to maintainer-clean-local if the HTML documentation should be part of the released tarball
clean-local:
	rm -rf cddlibman/

# Install the PDF, PS, and DVI documentation into PREFIX/share/doc/cddlib/
doc_DATA = cddlibman.pdf cddlibman.dvi cddlibman.ps
