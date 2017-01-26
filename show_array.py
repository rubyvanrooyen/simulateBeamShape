#! /usr/bin/python

import numpy
import os, string, sys
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Show any selected subarray
def showLayout(mkat, subarray=None, savegraph=False):
    [mkat_lat, mkat_lon, subarr_lat, subarr_lon] = buildArray(mkat.array_ref, mkat.array, subarray=mkat.__dict__[subarray])
    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(mkat_lon, mkat_lat, 'ko', alpha=0.3)
    plt.plot(subarr_lon, subarr_lat, 'ro', alpha=0.3)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.title(str(subarray))
    if savegraph: plt.savefig('SimpleSubArray.png')

def generateMap(mkat, subarray=None, antennalist=[], savegraph=False):
    if subarray is None and len(antennalist)<1: raise RuntimeError('No subarray specified. Provide config name or antenna list')
    if subarray is not None:
        if mkat.__dict__.has_key(subarray):
            [mkat_lat, mkat_lon, subarr_lat, subarr_lon] = buildArray(mkat.array_ref, mkat.array, subarray=mkat.__dict__[subarray], savegraph=savegraph)
        else: raise RuntimeError('Unknown array configuration %s' % subarray)
    else:
        [mkat_lat, mkat_lon, subarr_lat, subarr_lon] = buildArray(mkat.array_ref, mkat.array, subarray=antennalist, savegraph=savegraph)
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
    fig = plt.figure(figsize=(20,13), facecolor='white')
    ax = fig.add_subplot(111)
    # draw parallels.
    parallels = numpy.arange(-30.73,-30.69,.01)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    # draw meridians
    meridians = numpy.arange(21.41,21.48,.01)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    m.drawmapboundary(fill_color='white')
    m.scatter(ref_x,ref_y, 1000, marker='+',color='r', label='Array reference', alpha=0.3)
    m.scatter(mkat_x,mkat_y,3,marker='o',color='b', label='MeerKAT antennas', alpha=0.3)
    m.scatter(subarr_x,subarr_y,3, marker='o',color='k', label='%s antennas'%str(subarray))
    cntr=0
    for x,y in zip(subarr_x, subarr_y):
        plt.text(x,y,mkat.mkat[cntr], fontsize=6, ha='center', va='baseline', color='k')
        cntr+=1
    plt.title('Detail Map: SubArray Layout')
    plt.legend(loc=0, numpoints=1, ncol=1, prop={'size':10}, scatterpoints=1)
    if savegraph: plt.savefig('SubArrayLayout.png')

def buildArray(mkat_reference, mkat_antennas, subarray=[], savegraph=False):
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

# -fin-

