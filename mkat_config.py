#! /usr/bin/python
##
## Holds hardcoded and predefined parameters used by simulation scripts
##
###############

import string
import astropy
from astropy import units as u
from astropy.coordinates import Longitude, Latitude, EarthLocation

def read_array(mkatfile):
    fin = open(mkatfile,'r')
    data = fin.readlines()
    fin.close()
    # ignore first line with reference information, ref position repeated per
    # antenna location
    array = {}
    antennas={}
    for idx in range(1,len(data)):
        antenna = data[idx].strip().split(',')
        name = antenna[0].strip()
        if array.has_key(name): raise RuntimeError('Duplicate antenna name: Exiting')

        # Reference location
        ref_LAT = antenna[1].strip()
        ref_LON = antenna[2].strip()
        ref_ALT = antenna[3].strip()
        lon = Longitude(ref_LON.strip(), u.degree, wrap_angle=180*u.degree, copy=False)
        lat = Latitude(ref_LAT.strip(), u.degree, copy=False)
        height = u.Quantity(float(ref_ALT.strip()), u.m, copy=False)
        ref_location = EarthLocation(lat=lat.to(u.deg).value, lon=lon.to(u.deg).value, height=height.to(u.m).value)
        [x,y,z] = ref_location.to_geocentric()

        # Antenna location
        [delta_North, delta_East, delta_Up] = antenna[5].strip().split()
        array[name] = EarthLocation((x.value+float(delta_North))*u.m, (y.value+float(delta_East))*u.m, (z.value+float(delta_Up))*u.m)
        antennas[name] = [delta_North, delta_East, delta_Up]
    return [ref_location, array, antennas]

class Subarrays():
    def __init__(self, mkatfile):
        # Get array layout coordinates (64 antennas)
        [self.array_ref, self.array, self.antennas] = read_array(mkatfile)
        # AR array rollout plan 2015
        self.ar2 = ['m063', 'm062', 'm024', 'm025', 'm031', 'm034', 'm015', 'm014', 'm001', 'm003', 'm006', 'm010', 'm008', 'm007', 'm021', 'm022', 'm036', 'm017', 'm018', 'm020', 'm011', 'm012', 'm000', 'm002', 'm005', 'm042', 'm041', 'm040', 'm038', 'm037', 'm030', 'm028']
        # Standard subarray configurations
        self.ar1= [string.lower(ant) for ant in self.ar2[:16]]
        self.mkat=[string.lower(ant) for ant in self.array.keys()]
        self.mkatcore= ['m002', 'm000', 'm005', 'm006', 'm001', 'm003', 'm004', 'm018', 'm020', 'm017', 'm015', 'm029', 'm021', 'm007', 'm019', 'm009', 'm016', 'm011', 'm028', 'm012', 'm027', 'm034', 'm035', 'm014', 'm042', 'm022', 'm013', 'm036', 'm008', 'm031', 'm026', 'm047', 'm030', 'm010', 'm032', 'm041', 'm037', 'm023', 'm038', 'm043', 'm039', 'm040', 'm024', 'm025']

    def list_subs(self):
        subarrays=self.__dict__.keys()
        del subarrays[subarrays.index('array_ref')]
        del subarrays[subarrays.index('array')]
        del subarrays[subarrays.index('antennas')]
        print 'Known subarray configurations:\n\t%s' % ', '.join(subarrays)

    def check_sub(self, subname):
        if not self.__dict__.has_key(subname): return False
        else: return True

    def get_sub(self, subname):
        return self.__dict__[subname]

# -fin-
