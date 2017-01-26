#! /usr/bin/python
import os, sys
import numpy
import matplotlib.pyplot as plt

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
    import optparse
    usage='%prog [options]'
    parser = optparse.OptionParser(usage=usage, version="%prog 1.0")
    parser.add_option('--mkat',
                      action='store',
                      dest='mkat',
                      type=str,
                      default='config/mkat_antennas.txt',
                      help='MeerKAT antenna positions provided by Ludwig 2015')
    parser.add_option('--sub',
                      action='store',
                      dest='subarray',
                      type=str,
                      default='mkat',
                      help='Name of subarray as defined in config')
    parser.add_option("--map",
                      dest='genmap',
                      action="store_true",
                      default=False,
                      help="Generate detailed map of subarray layout")
    parser.add_option("--silent",
                      dest='silent',
                      action="store_true",
                      default=False,
                      help="Do not show any messages on stdout for scripting")
    parser.add_option('-o', "--output",
                      dest='savegraph',
                      action="store_true",
                      default=False,
                      help="Save graphs to PNG format")
    parser.add_option('-v', "--verbose",
                      dest='verbose',
                      action="store_true",
                      default=False,
                      help="Display intermittend results and all graphs")
    (opts, args) = parser.parse_args()

    from mkat_config import Subarrays
    mkat = Subarrays(opts.mkat)

    if opts.subarray=='?':
        mkat.list_subs()
        import sys
        sys.exit(0)
    if not mkat.check_sub(opts.subarray):
        print('\nUnknown subarray configuration \'%s\'' % opts.subarray)
        print('Available subarrays')
        mkat.list_subs()
        print
        raise RuntimeError('Unknown subarray configuration')

    import show_array
    if opts.genmap:
        show_array.generateMap(mkat, opts.subarray, savegraph=opts.savegraph)
    show_array.showLayout(mkat, opts.subarray, savegraph=opts.savegraph)

    # Get subarray layout coordinates
    filename = make_simms_array(mkat, opts.subarray)
    if not opts.silent: print 'SIMMS output file written to %s' % filename
    else: print filename

    if opts.verbose:
        try: plt.show()
        except: pass # nothing to show

# -fin-
