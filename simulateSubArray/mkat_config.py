#! /usr/bin/python
##
## Holds hardcoded and predefined parameters used by simulation scripts
##
###############

import string
import astropy
from astropy import units as u
from astropy.coordinates import Longitude, Latitude, EarthLocation


# MeerKAT antenna positions provided by Ludwig 2015
mkatfile = 'arrays/mkat_antennas.txt'
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

# AR array rollout plan 2015
ar2 = ['M063', 'M062', 'M024', 'M025', 'M031', 'M034', 'M015', 'M014', 'M001', 'M003', 'M006', 'M010', 'M008', 'M007', 'M021', 'M022', 'M036', 'M017', 'M018', 'M020', 'M011', 'M012', 'M000', 'M002', 'M005', 'M042', 'M041', 'M040', 'M038', 'M037', 'M030', 'M028']

class Subarrays():
    def __init__(self):
        # Get array layout coordinates (64 antennas)
        [self.array_ref, self.array, self.antennas] = read_array(mkatfile)
        self.ar1= [string.lower(ant) for ant in ar2[:16]]
        self.ar2= [string.lower(ant) for ant in ar2]
        self.m64full=[string.lower(ant) for ant in self.array.keys()]

        # Subarrays (fixed antenna selection)
        self.m64core= ['m002', 'm000', 'm005', 'm006', 'm001', 'm003', 'm004', 'm018', 'm020', 'm017', 'm015', 'm029', 'm021', 'm007', 'm019', 'm009', 'm016', 'm011', 'm028', 'm012', 'm027', 'm034', 'm035', 'm014', 'm042', 'm022', 'm013', 'm036', 'm008', 'm031', 'm026', 'm047', 'm030', 'm010', 'm032', 'm041', 'm037', 'm023', 'm038', 'm043', 'm039', 'm040', 'm024', 'm025']
        self.m32core= ['m039', 'm011', 'm018', 'm013', 'm036', 'm012', 'm038', 'm015', 'm040', 'm008', 'm006', 'm047', 'm009', 'm041', 'm005', 'm010', 'm021', 'm029', 'm043', 'm004', 'm023', 'm025', 'm022', 'm002', 'm030', 'm000', 'm001', 'm034', 'm031', 'm027', 'm032', 'm017']
        self.m32wide= ['m055', 'm054', 'm052', 'm056', 'm033', 'm050', 'm053', 'm044', 'm045', 'm063', 'm059', 'm062', 'm061', 'm058', 'm049', 'm046', 'm060', 'm057', 'm048', 'm051', 'm003', 'm020', 'm007', 'm019', 'm016', 'm028', 'm035', 'm014', 'm042', 'm026', 'm037', 'm024']
        self.m32gena0= ['m060', 'm048', 'm062', 'm058', 'm053', 'm046', 'm054', 'm045', 'm055', 'm050', 'm013', 'm039', 'm018', 'm037', 'm036', 'm011', 'm038', 'm005', 'm006', 'm040', 'm047', 'm010', 'm024', 'm029', 'm003', 'm043', 'm004', 'm032', 'm030', 'm027', 'm031', 'm017']
        self.m32genb0= ['m061', 'm057', 'm059', 'm049', 'm056', 'm063', 'm052', 'm044', 'm033', 'm051', 'm015', 'm020', 'm035', 'm016', 'm012', 'm014', 'm019', 'm042', 'm009', 'm007', 'm008', 'm041', 'm002', 'm021', 'm023', 'm022', 'm026', 'm001', 'm028', 'm000', 'm034', 'm025']
        self.m32gena1= ['m060', 'm048', 'm062', 'm058', 'm053', 'm046', 'm054', 'm045', 'm033', 'm051', 'm013', 'm039', 'm018', 'm037', 'm036', 'm011', 'm038', 'm005', 'm006', 'm040', 'm047', 'm010', 'm024', 'm029', 'm003', 'm043', 'm004', 'm032', 'm030', 'm027', 'm031', 'm017']
        self.m32genb1= ['m061', 'm057', 'm059', 'm049', 'm056', 'm063', 'm052', 'm044', 'm055', 'm050', 'm015', 'm020', 'm035', 'm016', 'm012', 'm014', 'm019', 'm042', 'm009', 'm007', 'm008', 'm041', 'm002', 'm021', 'm023', 'm022', 'm026', 'm001', 'm028', 'm000', 'm034', 'm025']

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
