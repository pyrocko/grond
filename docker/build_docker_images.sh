#!/bin/bash
docker build nest -t grond-nest
docker build docs -t grond-docs


if [ ! -d "fat-nest/grond-test-data" ] ; then
    if [ ! -d "../test/data" ] ; then
        echo "Make sure complete test data is in "../test/data", by running the tests."
        exit 1
    fi

    rsync -av "../test/data/" "fat-nest/grond-test-data/"
    rsync -av "../test/data/" "fat-util/grond-test-data/"

fi

docker build fat-nest -t grond-fat-nest
docker build fat-util -t grond-fat-util
