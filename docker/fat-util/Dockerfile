FROM pyrocko-util

# additional runtime requirements for gmt
RUN apt-get update
RUN apt-get install -y \
        gmt gmt-gshhg poppler-utils imagemagick

COPY grond-test-data /grond-test-data
