FROM grond-nest

# docs requirements
RUN apt-get update
RUN apt-get install -y \
    texlive-fonts-recommended texlive-latex-extra \
    texlive-latex-recommendedt exlive-generic-extra
RUN pip3 install sphinx sphinxcontrib-programoutput git+https://git.pyrocko.org/pyrocko/sphinx-sleekcat-theme.git
