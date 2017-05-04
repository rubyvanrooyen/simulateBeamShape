#! /usr/bin/python
##
## Show layout of selected subarray on MeerKAT array
##
###############

import numpy
import os, string, sys
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    import numpy as np
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.

    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")

    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2.  * np.arctan2(tu, cf)

    cu = 1.  / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1.  + np.sqrt(1.  + c2a * (1.  / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1.  - x
    c = (x * x / 4.  + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):

        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2.  * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4.  - 3.) * y * cz * d / 6.  + x) * d / 4.  - cz) * sy * d + tu

    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3.  * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1.  - c) * d * f + np.pi) % (2*np.pi)) - np.pi

    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)

    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi

    return (glon2, glat2, baz)

def equi(m, centerlon, centerlat, radius, *args, **kwargs):
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])

    #~ m.plot(X,Y,**kwargs)
    #Should work, but doesn't...
    X,Y = m(X,Y)
    plt.plot(X,Y,**kwargs)

# Show any selected subarray
def mkat_subarr(mkat, subarray=None, antennalist=[], savegraph=False):
    if subarray is None and len(antennalist)<1: raise RuntimeError('No subarray specified. Provide config name or antenna list')
    if subarray is not None:
        if mkat.__dict__.has_key(subarray):
            [mkat_lat, mkat_lon, subarr_lat, subarr_lon] = build_array(mkat.array_ref, mkat.array, subarray=mkat.__dict__[subarray], savegraph=savegraph)
        else: raise RuntimeError('Unknown array configuration %s' % subarray)
    else:
        [mkat_lat, mkat_lon, subarr_lat, subarr_lon] = build_array(mkat.array_ref, mkat.array, subarray=antennalist, savegraph=savegraph)
    # Earth projection
    m = Basemap(projection='merc',
                lat_0=mkat.array_ref.latitude.value,
                lon_0=mkat.array_ref.longitude.value,
                llcrnrlon=numpy.min(mkat_lon)-0.005,
                llcrnrlat=numpy.min(mkat_lat)-0.005,
                urcrnrlon=numpy.max(mkat_lon)+0.005,
                urcrnrlat=numpy.max(mkat_lat)+0.005)
    # set regular grid and map projected coordinates
    ref_x, ref_y=m(mkat.array_ref.longitude.value, mkat.array_ref.latitude.value)
    mkat_x,mkat_y=m(mkat_lon, mkat_lat)
    subarr_x,subarr_y=m(subarr_lon, subarr_lat)
    # Array layout
    fig = plt.figure(figsize=(20,13))
    ax = fig.add_subplot(111)
    # draw parallels.
    parallels = numpy.arange(-30.73,-30.69,.01)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    # draw meridians
    meridians = numpy.arange(21.41,21.48,.01)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    m.drawmapboundary(fill_color='white')
    m.scatter(ref_x,ref_y, 500, marker='+',color='r', label='Array reference')
    m.scatter(mkat_x,mkat_y,5,marker='o',color='c', label='MeerKAT antennas')
    m.scatter(subarr_x,subarr_y,10, marker='o',color='k', label='SubArray')
    equi(m, mkat.array_ref.longitude.value, mkat.array_ref.latitude.value, radius=0.5,lw=1., linestyle='--', label='0.5 deg')
    equi(m, mkat.array_ref.longitude.value, mkat.array_ref.latitude.value, radius=1,lw=1., linestyle='--', label='1 deg')
    equi(m, mkat.array_ref.longitude.value, mkat.array_ref.latitude.value, radius=2,lw=1., linestyle='--', label='2 deg')
    equi(m, mkat.array_ref.longitude.value, mkat.array_ref.latitude.value, radius=3,lw=1., linestyle='--', label='3 deg')
    plt.title('SubArray Layout')
    plt.legend(loc=0, numpoints=1)
    if savegraph: plt.savefig('SubArrayLayout.png')

# Show AR1 and AR2 layout
def mkat_ar2(mkat, savegraph=False):
    [mkat_lat, mkat_lon, subarr_lat, subarr_lon] = build_array(mkat.array_ref, mkat.array, subarray=mkat.ar2, savegraph=savegraph)
    # First array releases
    fig = plt.figure(figsize=(20,13))
    ax = fig.add_subplot(111)
    m = Basemap(projection='merc',
                llcrnrlon=numpy.min(subarr_lon)-0.005,
                llcrnrlat=numpy.min(subarr_lat)-0.005,
                urcrnrlon=numpy.max(subarr_lon)+0.005,
                urcrnrlat=numpy.max(subarr_lat)+0.005)
    subarr_x,subarr_y=m(subarr_lon, subarr_lat)
    # draw parallels.
    parallels = numpy.arange(numpy.min(subarr_lat)-0.005,numpy.max(subarr_lat)+0.005,.005)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    # draw meridians
    meridians = numpy.arange(numpy.min(subarr_lon)-0.005,numpy.max(subarr_lon)+0.005,.005)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    m.drawmapboundary(fill_color='white')
    m.scatter(subarr_x,subarr_y,5, marker='o',color='b', label='MeerKAT AR2')
    m.scatter(subarr_x[:16],subarr_y[:16],5, marker='o',color='r', label='MeerKAT AR1')
    cntr=0
    for x,y in zip(subarr_x[:16], subarr_y[:16]):
        plt.text(x,y,mkat.ar2[cntr], fontsize=10, fontweight='bold', ha='center', va='center', color='k')
        cntr+=1
    plt.title('MeerKAT AR2 Layout')
    plt.legend(loc=0, numpoints=1)
    if savegraph: plt.savefig('MeerKATARArraysLayout.png')

# Show the MeerKAT array layout
def mkat_rollout(mkat, savegraph=False):
    [mkat_lat, mkat_lon, subarr_lat, subarr_lon] = build_array(mkat.array_ref, mkat.array, subarray=mkat.ar2, savegraph=savegraph)
    # Earth projection
    m = Basemap(projection='merc',
                lat_0=mkat.array_ref.latitude.value,
                lon_0=mkat.array_ref.longitude.value,
                llcrnrlon=numpy.min(mkat_lon)-0.005,
                llcrnrlat=numpy.min(mkat_lat)-0.005,
                urcrnrlon=numpy.max(mkat_lon)+0.005,
                urcrnrlat=numpy.max(mkat_lat)+0.005)
    # set regular grid and map projected coordinates
    ref_x, ref_y=m(mkat.array_ref.longitude.value, mkat.array_ref.latitude.value)
    mkat_x,mkat_y=m(mkat_lon, mkat_lat)
    subarr_x,subarr_y=m(subarr_lon, subarr_lat)
    # Array layout
    fig = plt.figure(figsize=(20,13))
    ax = fig.add_subplot(111)
    # draw parallels.
    parallels = numpy.arange(-30.73,-30.69,.01)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    # draw meridians
    meridians = numpy.arange(21.41,21.48,.01)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    m.drawmapboundary(fill_color='white')
    m.scatter(ref_x,ref_y, 500, marker='+',color='r', label='Array reference')
    m.scatter(mkat_x,mkat_y,5,marker='o',color='k', label='MeerKAT antennas')
    m.scatter(subarr_x,subarr_y,5, marker='o',color='b', label='MeerKAT AR2')
    m.scatter(subarr_x[:16],subarr_y[:16],5, marker='o',color='r', label='MeerKAT AR1')
    equi(m, mkat.array_ref.longitude.value, mkat.array_ref.latitude.value, radius=0.5,lw=1., linestyle='--', label='0.5 deg')
    equi(m, mkat.array_ref.longitude.value, mkat.array_ref.latitude.value, radius=1,lw=1., linestyle='--', label='1 deg')
    equi(m, mkat.array_ref.longitude.value, mkat.array_ref.latitude.value, radius=2,lw=1., linestyle='--', label='2 deg')
    equi(m, mkat.array_ref.longitude.value, mkat.array_ref.latitude.value, radius=3,lw=1., linestyle='--', label='3 deg')
    plt.title('MeerKAT Array Layout')
    plt.legend(loc=0, numpoints=1)
    if savegraph: plt.savefig('MeerKATArrayLayout.png')

def build_array(mkat_reference, mkat_antennas, subarray=[], savegraph=False):
    mkat_latitude=[]
    mkat_longitude=[]
    for antenna in numpy.sort(mkat_antennas.keys()):
        mkat_latitude.append(mkat_antennas[antenna].latitude.value)
        mkat_longitude.append(mkat_antennas[antenna].longitude.value)

    # Get subarray
    subarr_latitude=[]
    subarr_longitude=[]
    for antenna in subarray:
        subarr_latitude.append(mkat_antennas[antenna].latitude.value)
        subarr_longitude.append(mkat_antennas[antenna].longitude.value)

    return [mkat_latitude, mkat_longitude, subarr_latitude, subarr_longitude]

if __name__ == '__main__':
# python build_array.py -o -v
# python build_array.py --ar2 -o -v
# python build_array.py --sub ar1 -v
    import optparse
    usage='\
\n  %prog [options] \
\nExamples \
\n\tShow MeerKAT array: %prog \
\n\tShow AR1 and AR2: %prog --ar2 -v\
\n\tShow custom array: %prog --sub <name> -o\
\n\tList all known subarrays: %prog --sub ?\
'
    parser = optparse.OptionParser(usage=usage, version="%prog 1.0")
    parser.add_option('--sub',
                      action='store',
                      dest='subarray',
                      type=str,
                      default=None,
                      help='Name of subarray as defined in config')
    parser.add_option("--mkat",
                      dest='mkat',
                      action="store_true",
                      default=False,
                      help="Show MeerKAT array layout with AR1 and AR2 highlighted")
    parser.add_option("--ar2",
                      dest='ar2',
                      action="store_true",
                      default=False,
                      help="Show AR1 and AR2")
    parser.add_option('-o', "--output",
                      dest='savegraph',
                      action="store_true",
                      default=False,
                      help="Save array layout graph to PNG format")
    parser.add_option('-v', "--verbose",
                      dest='verbose',
                      action="store_true",
                      default=False,
                      help="Display array layout graph")
    (opts, args) = parser.parse_args()
    if not opts.ar2 and opts.subarray is None: opts.mkat=True

    from mkat_config import Subarrays
    mkat = Subarrays()

    if opts.subarray=='?':
        mkat.list_subs()
        import sys
        sys.exit(0)

    # generate and display array
    if opts.ar2:
        mkat_ar2(mkat, savegraph=opts.savegraph)
    if opts.subarray is not None:
        mkat_subarr(mkat, subarray=opts.subarray, savegraph=opts.savegraph)
    if opts.mkat:
        mkat_rollout(mkat, savegraph=opts.savegraph)

    if opts.verbose: plt.show()

# -fin-

