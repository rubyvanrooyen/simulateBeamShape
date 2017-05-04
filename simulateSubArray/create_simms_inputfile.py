#! /usr/bin/python
import os, sys
import numpy

# A standard file should have the format: pos1 pos2 pos3* dish_diameter station mount
def make_simms_array(mkat, subarray=None):
    if subarray is None: subarray=mkat
    coord_file = '%s_antennas.enu' % os.path.splitext(os.path.basename(subarray))[0]
    fout = open(coord_file, 'w')
    fout.write('#E N U dish_diam station mount\n')
    for ant in numpy.sort(mkat.get_sub(subarray)):
        [delta_N, delta_E, delta_H] = mkat.antennas[ant]
        fout.write('%s %s %s 13.5 %s ALT-AZ\n' % (delta_N, delta_E, delta_H, ant))
    fout.close()

    return coord_file

if __name__ == '__main__':

    from mkat_config import Subarrays
    mkat = Subarrays()

    try:
        subarray = sys.argv[1]
    except:
        print ('Input filename required: %s <subarray>' % sys.argv[0])
        raise SystemExit(mkat.list_subs())
    if not mkat.check_sub(subarray): raise RuntimeError('Unknown subarray configuration %s' % subarray)

    # Get subarray layout coordinates
    filename = make_simms_array(mkat, subarray)
    print 'SIMMS output file written to %s' % filename

# -fin-
