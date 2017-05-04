

from simms import simms
for deg in [-90,-30,0]:
    msname = "test_%ddeg.ms_simms" % deg
    simms.create_empty_ms(msname=msname,
                          tel='MeerKAT',
                          pos='mkat_antennas.enu',
                          pos_type='ascii',
                          coords='enu',
                          ra='0h0m0s',
                          dec='%ddeg'%deg,
                          freq0=856.0e+06,
                          dfreq=1671872.0,
                          synthesis=12,
                          dtime=256,
                          stokes='XX XY YX YY',
                          date=["UTC,2014/01/15"],
                          setlimits=True,
                          elevation_limit=20,
                          nbands=1,
                          nchan=512,
                          auto_corr=True,
                         )

    break

print msname
from casac import *

im = casac.imager()
im.open(msname)
im.calcuvw()
im.plotuv(False)
im.done()




