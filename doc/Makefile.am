# The Makefile for doc/ is non-standard as we distribute a manual that is
# built from a .tex file (and not a .texi file.)
# We do not want to make latex & pdflatx & dvips requirements for building
# cddlib, so we ship the built manual as part of the distribution.

LATEX_SRC = $(srcdir)/cddlibman.tex $(srcdir)/cddlibman.bbl

# Building twice seems to be enough to get references right.
cddlibman.dvi: ${LATEX_SRC}
if LATEX
	$(latex) $<
	$(latex) $<
else
	@echo "LaTeX was not found. Please build from the released source tarball or install LaTeX and rerun ./configure."
	@false
endif

cddlibman.pdf: ${LATEX_SRC}
if PDF
	$(pdflatex) $<
	$(pdflatex) $<
else
	@echo "PdfLatex was not found. Please build from the released source tarball or install PdfLatex and rerun ./configure."
	false
endif

cddlibman.ps: cddlibman.dvi
if PS
	$(dvips) cddlibman.dvi
else
	@echo "dvips was not found. Please build from the released source tarball or install dvips and rerun ./configure."
	false
endif

# Wipe the HTML output directory as latex2html's perl calls otherwise try to
# find images.pl in there without putting cddlibman into @INC.
# We depend on cddlibman.dvi so we know that the references have been correctly
# resolved in the .aux file.
cddlibman/cddlibman.html: ${LATEX_SRC} html.sty cddlibman.dvi
if HTML
	rm -rf cddlibman
	$(latex2html) $<
else
	@echo "latex2html was not found. Please build from the released source tarball or install latex2html and rerun ./configure."
	false
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
