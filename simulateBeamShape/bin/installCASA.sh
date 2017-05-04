#! /bin/bash

echo "Installing new version of CASA"

casadir="/usr/local/casa"
while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -f|--installfile)
        casafile="$2"
        shift # past argument
        ;;
        -l|--location)
        casadir="$2"
        shift # past argument
        ;;
        *)
        # unknown option
        ;;
    esac
    shift # past argument or value
done
scriptname=`basename "$0"`
if [[ -z "$casafile" ]]
then
    echo "Usage: $scriptname -l </usr/local> -f <casa-release-4.4.0-el6.tar.gz>"
fi

# Find full path names
if [ ! -d "$casadir" ]
then
    # directory does not exist -- ask the user to make one
    echo "Directory $casadir does not exist"
    echo "Please create directory and ensure you have write permission"
fi
echo $casadir
pushd $casadir > /dev/null
    casadir=`pwd`
popd > /dev/null
echo $casadir
installdir=`dirname "$casadir"`
installdir="$installdir/bin"
echo $installdir

# Removing old binaries and cleaning up old casa
# check if directory is empty
if [ "$(ls -A $casadir)" ]
then
    echo "directory is not empty"
    currentversion=$(find $casadir -mindepth 1 -maxdepth 1 -type d)
    echo $currentversion
    # Remove old symbolic links
    pushd $installdir > /dev/null
        find $currentversion/ -maxdepth 1 -type f -perm /+x | while read ; do file=$(basename "$REPLY"); echo $file ; rm -f $file ; done
    popd > /dev/null
    # Package old casa version
    pushd $casadir > /dev/null
        tar -cvzf "$currentversion.tgz" $currentversion
        rm -rf $currentversion
    popd > /dev/null
fi

# Unpack new CASA release
cp $casafile $casadir
newversion=`basename "$casafile"`
pushd $casadir > /dev/null
    tar -xvzf $newversion
    rm $newversion
popd > /dev/null

# Verify old or newer version to get correct install
newdir=${newversion%.*.*}
version=$(awk -F- '{print $3}' <<< "$newdir")
version=${version%.*}
if [ $version \< 4.7 ]
then
    echo "old version"
    symlinkdir="$casadir/${newversion%.*.*}/"
else
    echo "new version"
    symlinkdir="$casadir/${newversion%.*.*}/bin/"
fi

# Create new symbolic links to install directory
pushd $installdir > /dev/null
    find $symlinkdir -maxdepth 1 -type f -perm /+x | while read ; do echo "$REPLY" ; ln -s "$REPLY" ; done
popd > /dev/null

# -fin-
