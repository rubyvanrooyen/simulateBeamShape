#! /usr/bin/python
##
## Simulates uv-tracks, synthesize beam natural weighting, Gaussian taper given
## the observatory latitude, the antennas positions, the hour angle range for the
## tracks and source declination
##
###############

import numpy
import numpy as np
import os, string, sys
import matplotlib.pyplot as plt

VERBOSE=False

# Standards and known constants
c = 2.998e8 # velocity of light in m/s
f = 1420e6 #Hz
l= c/f   #m

# The function evaluates baselines lenghts and angle between antennas
def baseline_angles(antennaPosition):
    #number of antennas
    na = len(antennaPosition)
    #number pair or baseline
    nbl = na*(na-1)/2

    length_angle = np.zeros((nbl, 2))
    k = 0
    for i in range(na):
        for j in range(i+1, na):
            length_angle[k,0] = np.sqrt((antennaPosition[i,0]-antennaPosition[j,0])**2 + (antennaPosition[i,1]-antennaPosition[j,1])**2)
            length_angle[k,1] = np.arctan2((antennaPosition[i,1]-antennaPosition[j,1]) , (antennaPosition[i,0]-antennaPosition[j,0]))
            k = k +1
    return length_angle

# The following function transform baseline to x,y,z coordinates
# lengthbaseline: Baseline length
# elevation: Elevation angle in radian
# azimuth: Azimuth angle in radian
# The result is a vector of (x,y,z) with unit the same as baseline length
# Interferometry and Synthesis in Radio Astronomy, Chapter 4, Equation 4.4
def baseline_to_xyz(lengthbaseline, elevation, azimuth, latitude):
    x = np.cos(latitude)*np.sin(elevation) - np.sin(latitude)*np.cos(elevation)*np.cos(azimuth)
    y = np.cos(elevation)*np.sin(azimuth)
    z = np.sin(latitude)*np.sin(elevation) + np.cos(latitude)*np.cos(elevation)*np.cos(azimuth)
    xyz = np.array([(x,y,z)])
    return lengthbaseline * xyz.T

# Transform x, y, z to u, v, w components
# ha: Source hour angle in radians
# dec: Source declination in radian
# Interferometry and Synthesis in Radio Astronomy, Chapter 4, Equation 4.1
def xyz_to_baseline(ha, dec):

    a1 = np.sin(ha)
    a2 = np.cos(ha)
    a3 = 0.

    b1 = -1*np.sin(dec)*np.cos(ha)
    b2 = np.sin(dec)*np.sin(ha)
    b3 = np.cos(dec)

    c1 = np.cos(dec)*np.cos(ha)
    c2 = -1*np.cos(dec)*np.sin(ha)
    c3 = np.sin(dec)

    return np.array([(a1,a2,a3),(b1,b2,b3),(c1,c2,c3)])

# Evaluate the track of a single antenna pair
# ntimeslots: is the number timeslots
# listha : is the hour angle range
# The function return a set of u,v,w components
def track_uv(listha, lengthbaseline, elevation, azimuth, latitude, dec, ntimeslots):
    UVW = np.zeros((ntimeslots, 3), dtype=float)
    for i in range(ntimeslots):
        UVW[i, :] = np.dot(xyz_to_baseline(listha[i], dec),baseline_to_xyz(lengthbaseline, azimuth, elevation, latitude)).T
    return UVW

def uvMask(ha, lengthbaseline, elevation, azimuth, latitude, dec, ntimeslots, maxsize, uvscaling):
    maskmat = np.zeros((maxsize,maxsize))
    uvw = track_uv(ha, lengthbaseline, elevation, azimuth, latitude, dec, ntimeslots);
    sctrl = maxsize/2 + 1;
    #print uvw.shape
    for i in range(ntimeslots):
        maskmat[sctrl+int(np.ceil(uvw[i,0]*uvscaling)),sctrl+int(np.ceil(uvw[i,1]*uvscaling))] = 1.
        maskmat[sctrl-int(np.ceil(uvw[i,0]*uvscaling)),sctrl-int(np.ceil(uvw[i,1]*uvscaling))] = 1.
    return maskmat

# Gaussian Kernel for tapering given FWHM
# this can be done in a more beautiful way with numpy arrays
def Gauss (npix, FWHM):
    n = npix/2.
    scale = -1/(2*FWHM**2)
    G = np.zeros((npix,npix))
    for i in range(npix):
        for j in range(npix):
            G[i][j] = np.exp(scale*((i-n)**2 + (j-n)**2))
    return G


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
    P = P/l # baseline
    mx = numpy.abs(P).max()

    plt.figure(figsize=(20,13))
    plt.clf()
    plt.hold(True)
    plt.plot(P[:,0], P[:,1], 'o', markersize=2)
    plt.xlim(-mx, mx)
    plt.ylim(-mx, mx+5)
    plt.axis('equal')
    plt.xlabel('E')
    plt.ylabel('N')
    plt.title('Antenna positions')
    plt.gca().invert_xaxis()
    plt.hold(False)

    # number of time slots
    ntimeslots = 300;
    # hour angle range in hours
    ha = np.linspace(-4, 4, ntimeslots)*np.pi/12;

    # source, obs postions
    #latitute
    latitute = (-30. + 42./60. + 48./3600.)*np.pi/180.

    # declination convert in radian
    declinations=[-60., -30., -10., 0., 10., 30.] # degrees
    for decidx in range(len(declinations)):
        dec_deg=declinations[decidx]
        dec = dec_deg*np.pi/180.

        # matrix to store baseline lengths and azimuth angles
        B = baseline_angles(P)
        #number of antennas
        na = np.size(P, 0)
        #number pair or baseline
        nbl = na*(na-1)/2

##Show UV tracks
        # Plot the uv-Coverage
        plt.figure(figsize=(20,13))
        plt.hold(True)
        for i in range(nbl):
            uv = track_uv(ha, B[i, 0], 0., B[i, 1], latitute, dec, ntimeslots)
            plt.plot(uv[:,0], uv[:,1], 'b.', markersize=1)
            plt.plot(-uv[:,0], -uv[:,1], 'r.', markersize=1)
        plt.xlabel('u')
        plt.ylabel('v')
        plt.title('Declination %f deg: UV coverage'%dec_deg)
        mb = 5*np.sqrt((uv**2).sum(1)).max();
        plt.xlim(-mb,mb);
        plt.ylim(-mb,mb);
        plt.axis('equal')
        plt.gca().invert_xaxis()
        plt.hold(False)
        plt.savefig('%s_uv_%d.png'%(subarray,int(dec_deg)))

##Show UV mask and synthesize beam, natural weighting without a taper
        # number of pixels
        npix = 2**10 # should be a power of two
        # kernel matrix
        mask = np.zeros((npix, npix))
        # uv-grid scale factor to fit tracks into matrix
        uvscale = npix/2/mb * 0.95 * 0.5
        for i in range(nbl):
            mask = mask + uvMask(ha, B[i, 0], 0., B[i, 1], latitute, dec, ntimeslots, npix, uvscale)
        psf = np.fft.ifftshift(np.fft.ifft2(mask.T)).real

        # plot the mask
        plt.figure(figsize=(20,13))
        plt.subplot(121)
        plt.imshow(mask.T)
        plt.xlabel('u')
        plt.ylabel('v')
        plt.title('uv coverage, "natural weighting"')
        plt.subplot(122)
        plt.imshow(psf[len(psf)/2-2**8:len(psf)/2+2**8, len(psf)/2-2**8:len(psf)/2+2**8], interpolation='nearest', origin='lower', vmin=psf.min()/1000., vmax=psf.max())
        plt.xlabel('u')
        plt.ylabel('v')
        plt.title('psf, "natural weighting"')

##Show UV mask and synthesize beam, natural weighting with Gaussian taper
        gauss_kernel = Gauss (npix, npix/20.)
        psf = np.fft.ifftshift(np.fft.ifft2(gauss_kernel*(mask.T))).real

        # plot the mask
        plt.figure(figsize=(20,13))
        plt.subplot(121)
        plt.imshow(gauss_kernel*(mask.T))
        plt.xlabel('u')
        plt.ylabel('v')
        plt.title('uv coverage, "tapered natural weighting"')
        plt.subplot(122)
        plt.imshow(psf[len(psf)/2-2**8:len(psf)/2+2**8, len(psf)/2-2**8:len(psf)/2+2**8], interpolation='nearest', origin='lower', vmin=psf.min()/1000., vmax=psf.max())
        plt.xlabel('u')
        plt.ylabel('v')
        plt.title('psf, "tapered Gaussian weighting"')
        plt.savefig('%s_psf_%d.png'%(subarray,int(dec_deg)))

        # plot the psf
        plt.figure(figsize=(20,13))
        plt.imshow(psf[len(psf)/2-2**8:len(psf)/2+2**8, len(psf)/2-2**8:len(psf)/2+2**8], interpolation='nearest', origin='lower', vmin=psf.min()/1000., vmax=psf.max())
        plt.xlabel('u')
        plt.ylabel('v')
        plt.title('Declination %f deg: PSF, "tapered Gaussian weighting"'%dec_deg)
        plt.savefig('%s_psf_%d.png'%(subarray,int(dec_deg)))

    if VERBOSE: plt.show()





# -fin-

