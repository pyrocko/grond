# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = grond
SOURCEDIR     = source
BUILDDIR      = build

GRONDCOMMANDS = scenario init events check go continue forward harvest plot movie \
   export report diff qc-polarization upgrade-config version 

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

source/cli/scenario.rst:
	for cmd in $(GRONDCOMMANDS) ; do \
		grond $$cmd --docs > source/cli/$${cmd}.rst ; \
	done


clean:
	for cmd in $(GRONDCOMMANDS) ; do \
		rm source/cli/$${cmd}.rst ; \
	done
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile source/cli/scenario.rst
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
