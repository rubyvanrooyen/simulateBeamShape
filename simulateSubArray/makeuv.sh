#! /bin/bash

for subarr in ar1 ar2 m32core m32wide m32gena1 m32gena0 m32genb0 m32genb1 m64core m64full
do
  echo $subarr
  python uv_coverage.py $subarr
done

# -fin-
