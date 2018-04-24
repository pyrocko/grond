#!/bin/bash
set -e
VERSION=v`python3 -c "import grond; print(grond.__version__);"`

if [ ! -f maintenance/deploy-docs.sh ] ; then
    echo "must be run from pyrocko's toplevel directory"
    exit 1
fi

cd doc
rm -rf build/$VERSION
make clean; make html $1
cp -r build/html build/$VERSION

read -r -p "Are your sure to update live docs at http://pyrocko.org/docs/grond/$VERSION [y/N]?" resp
case $resp in
    [yY][eE][sS]|[yY] )
        rsync -av build/$VERSION pyrocko@hive:/var/www/pyrocko.org/docs/grond/;
        ;;
    * ) ;;
esac

read -r -p "Do you want to link 'current' to the just uploaded version $VERSION [y/N]?" resp
case $resp in
    [yY][eE][sS]|[yY] )
        echo "Linking docs/grond/$VERSION to docs/grond/current";
        ssh pyrocko@hive "rm -rf /var/www/pyrocko.org/docs/grond/current; ln -s /var/www/pyrocko.org/docs/grond/$VERSION /var/www/pyrocko.org/docs/grond/current";
        ;;
    * ) ;;
esac

cd ..
