FROM pyrocko

WORKDIR /src
RUN pip3 install utm
RUN git clone https://github.com/pyrocko/kite.git && cd kite \
    && python3 setup.py install
