# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

DOXYDOC_OUTDIR=$(BUILDDIR)/../../build/doc/html/sphinx


# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	$(shell [[ -d "$(DOXYDOC_OUTDIR)" ]] && rm -Rf "$(DOXYDOC_OUTDIR)")
	$(shell cp -Rf $(BUILDDIR)/html $(DOXYDOC_OUTDIR))
	$(shell [[ -d "$(DOXYDOC_OUTDIR)/_sources" ]] && rm -Rf "$(DOXYDOC_OUTDIR)/_sources")

