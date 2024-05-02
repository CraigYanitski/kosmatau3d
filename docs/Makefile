# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = _source
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
# It also clears all generated files and 
# makes Sphinx think it is run on READTHEDOCS
%: Makefile
	rm -rf _build
	rm -rf _source/_autosummary
	rm -rf _source/generated
	rm -f _source/kosmatau3d.*
	jupytext --to notebook notebooks/*.py
	if [ -z "$$READTHEDOCS" ]; then export READTHEDOCS="temp"; else echo "READTHEDOCS"; fi
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	unset READTHEDOCS
