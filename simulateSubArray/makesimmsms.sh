#! /bin/bash

for subarr in ar1 ar2 m32core m32wide m32gena1 m32gena0 m32genb0 m32genb1 m64core m64full
do
  echo $subarr
  python make_cfg.py "configurations/"$subarr"_ANTENNAS" --cfg $subarr".cfg" --msname $subarr".MS"
  makems $subarr".cfg"
  mv $subarr".MS_p0/" $subarr".MS/"
done

# -fin-
