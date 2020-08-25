FROM grond-nest

# docs requirements
RUN apt-get install -y python3-sphinx \
    texlive-fonts-recommended texlive-latex-extra \
    texlive-latex-recommended texlive-generic-extra \
    python3-sphinxcontrib.programoutput

RUN pip3 install git+https://git.pyrocko.org/pyrocko/sphinx-sleekcat-theme.git
