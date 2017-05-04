#! /usr/bin/python
##
## Calculate longest and shortest baseline
##
###############

import numpy
import numpy as np
import os, string, sys
import matplotlib.pyplot as plt

if __name__ == '__main__':

    from mkat_config import Subarrays
    mkat = Subarrays()

    try:
        subarray = sys.argv[1]
    except:
        print ('Input filename required: %s <subarray>' % sys.argv[0])
        raise SystemExit(mkat.list_subs())
    if not mkat.check_sub(subarray): raise RuntimeError('Unknown subarray configuration %s' % subarray)

    delta_pos=[]
    for ant in numpy.sort(mkat.get_sub(subarray)):
        [delta_N, delta_E, delta_H] = mkat.antennas[ant]
        delta_pos.append([delta_N, delta_E])
    P=numpy.array(delta_pos,dtype='float')

    from uv_coverage import baseline_angles
    B = baseline_angles(P)
    bl_idx=numpy.argsort(B[:,0])
    print 'Subarray ', subarray
    print 'Shortest baseline %.2f m' % B[bl_idx[0],0]
    print 'Longest baseline %.2f km' % (B[bl_idx[-1],0]/1000.)

# -fin-

