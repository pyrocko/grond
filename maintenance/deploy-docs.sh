#!/bin/bash
set -e
VERSION=v`python3 -c "import grond; print(grond.__version__);"`

if [ ! -f maintenance/deploy-docs.sh ] ; then
    echo "must be run from grond's toplevel directory"
    exit 1
fi

cd docs
rm -rf build/$VERSION
make clean; make html $1
cp -r build/html build/$VERSION

read -r -p "Are your sure to update live docs at http://pyrocko.org/grond/docs/$VERSION/ [y/N]?" resp
case $resp in
    [yY][eE][sS]|[yY] )
        rsync -av build/$VERSION/ pyrocko@hive:/var/www/pyrocko.org/grond/docs/$VERSION/;
        ;;
    * ) ;;
esac

read -r -p "Do you want to link 'current' to the just uploaded version $VERSION [y/N]?" resp
case $resp in
    [yY][eE][sS]|[yY] )
        echo "Linking grond/docs/$VERSION to grond/docs/current";
        ssh pyrocko@hive "rm -f /var/www/pyrocko.org/grond/docs/current; ln -s /var/www/pyrocko.org/grond/docs/$VERSION /var/www/pyrocko.org/grond/docs/current";
        ;;
    * ) ;;
esac

cd ..
