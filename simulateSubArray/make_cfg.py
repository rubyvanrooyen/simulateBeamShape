#! /usr/bin/python

from optparse import OptionParser

def make_cfg(ANTENNAS, cfgname='makems.cfg', msname='simple.ms'):
    print msname
    fout = open(cfgname, 'w')
    fout.write('AntennaTableName=%s\n' % ANTENNAS) # list of antenna positions (CASA table)
    fout.write('NTimes=300\n') # number of integrations
    # fout.write('StepTime=60.\n') # integration time in seconds
    fout.write('StepTime=86400.\n') # integration time in seconds
    fout.write('StartFreq=856.0e+06\n')
    fout.write('StepFreq=1671872.0\n')
    fout.write('NFrequencies=512\n')
    fout.write('NBands=1\n')
    fout.write('Declination=-30:00:00\n')
    fout.write('RightAscension=00:00:00\n')
    fout.write('StartTime=2014/01/01/15:00\n')
    fout.write('MSName=%s\n' % msname)
    fout.write('MSDesPath=.\n') # Where to save MS
    fout.write('NParts=1\n')
    fout.write('WriteImagingColumns=True\n') # include imaging columns; {MODEL,CORRECTED}_DATA
    fout.write('WriteAutoCorr=True\n') # Include autocorrelations
    fout.close()

if __name__ == '__main__':
    usage = "\
python %prog [options] <fullpath/ANTENNAS> \
\n \
\nExample: \
\n\t python %prog -p ../output/ -a ../output/w2332_.kat7 -r 0.1 --sumss \
"
    parser = OptionParser(usage=usage, version="%prog 1.0")
    parser.add_option('--cfg',
                      action='store',
                      dest='cfgname',
                      type=str,
                      default=None,
                      help='Config output name')
    parser.add_option('--msname',
                      action='store',
                      dest='msname',
                      type=str,
                      default=None,
                      help='Name for output CASA MS')

    (opts, args) = parser.parse_args()
    if len(args)<1: raise SystemExit(parser.print_help())

    if opts.cfgname is None and opts.msname is None:
        make_cfg(args[0])
    elif opts.cfgname is not None and opts.msname is None:
        make_cfg(args[0], cfgname=opts.cfgname)
    elif opts.cfgname is None and opts.msname is not None:
        make_cfg(args[0], msname=opts.msname)
    else:
        make_cfg(args[0], cfgname=opts.cfgname, msname=opts.msname)

# -fin-
