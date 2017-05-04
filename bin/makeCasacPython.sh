#! /bin/bash -x

casa="/home/ruby/ska_development/simulateBeamShape/venv/local/casa/casa-release-4.5.3-el6"
python="/home/ruby/ska_development/simulateBeamShape/venv/local/lib/python2.7/site-packages"

pushd $casa/lib > /dev/null
ls
    # copy basic Python files
    ls $python
    cp -a python2.7/casac.py python2.7/__casac__ $python
    ls $python

    # copy dependent libraries, with moderate sophistication
    for f in lib*.so*
    do
        echo $f
        if [ -h $f ]
        then
            cp -a $f $python/__casac__ # copy symlinks as links
        else
            case $f in
                *_debug.so) ;; # skip -- actually text files
                libgomp.*)
                    # somehow patchelf fries this particular file
                    cp -a $f $python/__casac__ ;;
                *)
                    cp -a $f $python/__casac__
                    patchelf --set-rpath '$ORIGIN' $python/__casac__/$f ;;
            esac
        fi
    done
popd > /dev/null

# patch rpaths of Python module binary files
pushd $python/__casac__ > /dev/null
    for f in _*.so
    do
        patchelf --set-rpath '$ORIGIN' $f
    done
popd > /dev/null

# -fin-
