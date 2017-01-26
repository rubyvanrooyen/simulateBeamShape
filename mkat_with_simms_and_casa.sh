#! /bin/bash
now=$(date)
echo
echo
echo "Running SIMMS beamshape simulation: $now"
echo

while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -s|--subarr)
        subarr="$2"
        shift # past argument
        ;;
        *)
                # unknown option
        ;;
    esac
    shift # past argument or value
done

script=`basename "$0"`

if [[ -z "$subarr" ]]
then
    echo "Usage: $script -s <subarray>"
    exit 1
fi

# Remove older files and generated measurement sets"
make clobber

echo "Generating simulation files for $subarr"
# Create an input file for simms
ENUfile=$(python mkat_simms_setup.py --sub mkatcore --silent)
# Creating a CASA antenna table
simms -T MeerKAT -n $subarr".MS" -t ascii -cs enu -l mkat_antennas -dec -30d42m48s -ra 0h0m0s -st 1 -dt 60 $ENUfile
# Generate a configuration file
python make_cfg.py $subarr".MS/ANTENNA" --cfg $subarr".cfg" --msname $subarr".MS"
# Build CASA measurement set
makems $subarr".cfg"
mv $subarr".MS_p0/" $subarr".MS/"

echo
echo
echo "Run CASA to simulate beam"
echo "casa --log2term"
echo "clean(vis='$subarr".MS"',niter=0, imagename='$subarr')"
echo "imview(raster={'file': '$subarr".psf"', 'colormap': 'Misc. 1: Isophotes', 'scaling': -1}, out='$subarr"_psf.png"')"
echo

# -fin-
