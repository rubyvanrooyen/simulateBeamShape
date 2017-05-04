#! /bin/bash

for subarr in ar1 ar2 m32core m32wide m32gena1 m32gena0 m32genb0 m32genb1 m64core m64full
do
  echo $subarr
  python create_simms_inputfile.py $subarr
  simms -T MeerKAT -n $subarr".MS" -t ascii -cs enu -l $subarr"_antennas" -dec -30d42m48s -ra 0h0m0s -st 1 -dt 60 $subarr"_antennas.enu"
  mv $subarr".MS/ANTENNA/" "configurations/"$subarr"_ANTENNAS"
  rm -rf $subarr".MS"
done

# -fin-
