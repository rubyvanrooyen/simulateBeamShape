#! /usr/bin/python
import os, sys
import numpy
import matplotlib.pyplot as plt
import random

if __name__ == '__main__':

    from mkat_config import Subarrays
    mkat = Subarrays()
    name=[]
    longitude=[]
    latitude=[]
    radius=[]
    for antenna in mkat.array.keys():
        name.append(antenna)
        longitude.append(mkat.array[antenna].longitude.value)
        latitude.append(mkat.array[antenna].latitude.value)
        radius.append(numpy.sqrt((longitude[-1]-mkat.array_ref.longitude.value)**2+(latitude[-1]-mkat.array_ref.latitude.value)**2))
    name = numpy.array(name)
    radius = numpy.array(radius)
    latitude = numpy.array(latitude)
    longitude = numpy.array(longitude)

    # divide the array in sections to extract antennas for the subarrays from
    sortedradius = numpy.argsort(radius)
    core = sortedradius[:-19]
    bwide = sortedradius[-19:]
    bmax = sortedradius[-10:]
    # M051: name[core][-3] == add to wide and subtract form core to get equal number of antennas
    antidx = numpy.where(name ==name[core][-3])[0][0]
    core = numpy.delete(core, numpy.where(core==antidx)[0][0])
    bwide = numpy.hstack((bwide, antidx))

    # divide the array in sections to ensure antenna selection across the array
    nw = numpy.nonzero(numpy.logical_and((longitude-mkat.array_ref.longitude.value)<0 , (latitude-mkat.array_ref.latitude.value)>0))[0]
    cnw = numpy.nonzero(numpy.logical_and((longitude[core]-mkat.array_ref.longitude.value)<0 , (latitude[core]-mkat.array_ref.latitude.value)>0))[0]
    sw = numpy.nonzero(numpy.logical_and((longitude-mkat.array_ref.longitude.value)<0 , (latitude-mkat.array_ref.latitude.value)<0))[0]
    csw = numpy.nonzero(numpy.logical_and((longitude[core]-mkat.array_ref.longitude.value)<0 , (latitude[core]-mkat.array_ref.latitude.value)<0))[0]
    se = numpy.nonzero(numpy.logical_and((longitude-mkat.array_ref.longitude.value)>0 , (latitude-mkat.array_ref.latitude.value)<0))[0]
    cse = numpy.nonzero(numpy.logical_and((longitude[core]-mkat.array_ref.longitude.value)>0 , (latitude[core]-mkat.array_ref.latitude.value)<0))[0]
    ne = numpy.nonzero(numpy.logical_and((longitude-mkat.array_ref.longitude.value)>0 , (latitude-mkat.array_ref.latitude.value)>0))[0]
    cne = numpy.nonzero(numpy.logical_and((longitude[core]-mkat.array_ref.longitude.value)>0 , (latitude[core]-mkat.array_ref.latitude.value)>0))[0]

    # core array of 32 antennas
    core32=numpy.hstack((core[random.sample(cnw,8)], core[random.sample(csw,8)], core[random.sample(cse,8)], core[random.sample(cne,8)]))
    # all long baselines with 12 core antennas
    rem = numpy.delete(core, [numpy.where(core==idx)[0][0] for idx in core32])
    wide32=numpy.hstack((bwide,rem))

    # Divide the longer baselines to create 2 imaging subarrays
    baselines={}
    for idx1 in range(len(bwide)):
        antenna1 = name[bwide[idx1]]
        longitude1 = mkat.array[antenna1].longitude.value
        latitude1 = mkat.array[antenna1].latitude.value
        for idx2 in range(idx1+1,len(bwide)):
            antenna2 = name[bwide[idx2]]
            longitude2 = mkat.array[antenna2].longitude.value
            latitude2 = mkat.array[antenna2].latitude.value
            baselines[(numpy.sqrt((longitude1-longitude2)**2+(latitude1-latitude2)**2))]=[antenna1,antenna2]
    used=[]
    blen = numpy.sort(baselines.keys())[::-1]
    longA=[]
    longB=[]
    for idx in range(len(blen)):
        bl=blen[idx]
        # if not used assign to an array
        if (baselines[bl][0] in used) or (baselines[bl][1] in used): continue
        if len(longA) == len(longB):
            longA.extend(baselines[bl])
        elif len(longA) < len(longB):
            longA.extend(baselines[bl])
        else:
            longB.extend(baselines[bl])
        # remove from avialable
        used.extend(baselines[bl])

    # Randomly select the shorter baselines to fill 32 antenna arrays
    shortB=[]
    shortA=[]
    rem=[]
    nwidx=random.sample(cnw,len(cnw)/2)
    shortB.extend(nwidx)
    tmp=numpy.delete(cnw, [numpy.where(cnw==idx)[0][0] for idx in nwidx])
    shortA.extend(random.sample(tmp,len(cnw)/2))

    swidx=random.sample(csw,len(csw)/2)
    shortA.extend(swidx)
    tmp=numpy.delete(csw, [numpy.where(csw==idx)[0][0] for idx in swidx])
    shortB.extend(random.sample(tmp,len(csw)/2))

    seidx=random.sample(cse,len(cse)/2)
    shortB.extend(seidx)
    tmp=numpy.delete(cse, [numpy.where(cse==idx)[0][0] for idx in seidx])
    shortA.extend(random.sample(tmp,len(cse)/2))

    neidx=random.sample(cne,len(cne)/2)
    shortA.extend(neidx)
    tmp=numpy.delete(cne, [numpy.where(cne==idx)[0][0] for idx in neidx])
    shortB.extend(random.sample(tmp,len(cne)/2))

    shortA=core[shortA]
    shortB=core[shortB]
    rem = numpy.delete(core, [numpy.where(core==idx)[0][0] for idx in numpy.hstack((shortA,shortB))])
    shortA = numpy.hstack((shortA, rem[0]))
    shortB = numpy.hstack((shortB, rem[1]))

    # Display results for verification
    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(longitude, latitude, 'ko', alpha=0.3)
    plt.plot(longitude[nw], latitude[nw], 'ro', alpha=0.3)
    plt.plot(longitude[core][cnw], latitude[core][cnw], 'ro', alpha=0.5)
    plt.plot(longitude[sw], latitude[sw], 'bo', alpha=0.3)
    plt.plot(longitude[core][csw], latitude[core][csw], 'bo', alpha=0.5)
    plt.plot(longitude[se], latitude[se], 'go', alpha=0.3)
    plt.plot(longitude[core][cse], latitude[core][cse], 'go', alpha=0.5)
    plt.plot(longitude[ne], latitude[ne], 'yo', alpha=0.3)
    plt.plot(longitude[core][cne], latitude[core][cne], 'yo', alpha=0.5)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.title('Array Quadrants')

    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(longitude, latitude, 'bo', alpha=0.3)
    plt.plot(longitude[core], latitude[core], 'ro', alpha=0.3)
    cntr=0
    for x,y in zip(longitude, latitude):
        if name[cntr] in name[bwide]:
            plt.text(x,y,name[cntr], fontsize=10, fontweight='bold', ha='center', va='center', color='k')
        cntr+=1
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)

    # Build subarrays
    print "Full Array"
    print name
    print "Full Core Array"
    print name[core]
    print "Core of 32 antennas"
    print name[core32]
    print "All long baselines and 12 core antennas"
    print name[wide32]

    longA1=['m060', 'm048', 'm062', 'm058', 'm053', 'm046', 'm054', 'm045', 'm033', 'm051']
    longB1=['m061', 'm057', 'm059', 'm049', 'm056', 'm063', 'm052', 'm044', 'm055', 'm050']
    arrA0=numpy.hstack((longA, name[shortA]))
    print 'General 32 antenna array A0:'
    print arrA0
    arrB0=numpy.hstack((longB, name[shortB]))
    print 'General 32 antenna array B0:'
    print arrB0
    arrA1=numpy.hstack((longA1, name[shortA]))
    print 'General 32 antenna array A1:'
    print arrA1
    arrB1=numpy.hstack((longB1, name[shortB]))
    print 'General 32 antenna array B1:'
    print arrB1

    # Display all the Arrays
    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(longitude, latitude, 'ro', alpha=0.3)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.title('M64Full')
    plt.savefig('M64Full.png')

    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(longitude, latitude, 'ko', alpha=0.3)
    plt.plot(longitude[core], latitude[core], 'ro', alpha=0.3)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.title('M64Core')
    plt.savefig('M64Core.png')

    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(longitude, latitude, 'ko', alpha=0.3)
    plt.plot(longitude[core32], latitude[core32], 'ro', alpha=0.3)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.title('M32Core')
    plt.savefig('M32Core.png')

    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(longitude, latitude, 'ko', alpha=0.3)
    plt.plot(longitude[wide32], latitude[wide32], 'ro', alpha=0.3)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.title('M32Wide')
    plt.savefig('M32Wide.png')

    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(longitude, latitude, 'ko', alpha=0.3)
    plt.plot(longitude[core32], latitude[core32], 'ro', alpha=0.3)
    plt.plot(longitude[wide32], latitude[wide32], 'bo', alpha=0.3)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.title('M32Core and M32Wide')
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.savefig('M32CoreandM32Wide.png')

    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(longitude, latitude, 'ko', alpha=0.3)
    plt.plot(longitude[[numpy.where(name==ant)[0][0] for ant in arrA0]], latitude[[numpy.where(name==ant)[0][0] for ant in arrA0]], 'ro', alpha=0.3)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.title('M32Gen1 v1')
    plt.savefig('M32Gen1v1.png')

    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(longitude, latitude, 'ko', alpha=0.3)
    plt.plot(longitude[[numpy.where(name==ant)[0][0] for ant in arrB0]], latitude[[numpy.where(name==ant)[0][0] for ant in arrB0]], 'ro', alpha=0.3)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.title('M32Gen2 v1')
    plt.savefig('M32Gen2v1.png')

    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(longitude, latitude, 'ko', alpha=0.3)
    plt.plot(longitude[[numpy.where(name==ant)[0][0] for ant in arrA0]], latitude[[numpy.where(name==ant)[0][0] for ant in arrA0]], 'ro', alpha=0.3)
    plt.plot(longitude[[numpy.where(name==ant)[0][0] for ant in arrB0]], latitude[[numpy.where(name==ant)[0][0] for ant in arrB0]], 'bo', alpha=0.3)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.title('M32Gen1 v1 and M32Gen2 v1')
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.savefig('M32Gen1v1andM32Gen2v1.png')

    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(longitude, latitude, 'ko', alpha=0.3)
    plt.plot(longitude[[numpy.where(name==ant)[0][0] for ant in arrA1]], latitude[[numpy.where(name==ant)[0][0] for ant in arrA1]], 'ro', alpha=0.3)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.title('M32Gen1 v2')
    plt.savefig('M32Gen1v2.png')

    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(longitude, latitude, 'ko', alpha=0.3)
    plt.plot(longitude[[numpy.where(name==ant)[0][0] for ant in arrB1]], latitude[[numpy.where(name==ant)[0][0] for ant in arrB1]], 'ro', alpha=0.3)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.title('M32Gen2 v2')
    plt.savefig('M32Gen2v2.png')

    plt.figure(facecolor='white')
    ax1 = plt.axes(frameon=False)
    plt.plot(longitude, latitude, 'ko', alpha=0.3)
    plt.plot(longitude[[numpy.where(name==ant)[0][0] for ant in arrA1]], latitude[[numpy.where(name==ant)[0][0] for ant in arrA1]], 'ro', alpha=0.3)
    plt.plot(longitude[[numpy.where(name==ant)[0][0] for ant in arrB1]], latitude[[numpy.where(name==ant)[0][0] for ant in arrB1]], 'bo', alpha=0.3)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.title('M32Gen1 v2 and M32Gen2 v2')
    ax1.axes.get_yaxis().set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    plt.savefig('M32Gen1v2andM32Gen2v2.png')

    plt.show()


# -fin-

