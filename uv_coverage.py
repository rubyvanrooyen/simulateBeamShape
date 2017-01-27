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

# Display simulated beamshape
def plot_mask(mask, psf, filename= None, comment=''):
    plt.figure(figsize=(20,13), facecolor='white')
    plt.subplot(121)
    [nr, nc] = numpy.shape(mask)
    plt.imshow(mask)
    plt.xlim(0.5*nr-0.25*nr, 0.5*nr+0.25*nr)
    plt.ylim(0.5*nc-0.25*nc, 0.5*nc+0.25*nc)
    plt.xlabel('u')
    plt.ylabel('v')
    plt.title('uv coverage, "%s"'%comment)
    plt.subplot(122)
    [nr, nc] = numpy.shape(psf)
    plt.imshow(psf,
            interpolation='gaussian',
            origin='lower',
            vmin=psf.min()/1000.,
            vmax=psf.max())
    plt.xlim(0.5*nr-0.25*nr, 0.5*nr+0.25*nr)
    plt.ylim(0.5*nc-0.25*nc, 0.5*nc+0.25*nc)
    plt.xlabel('u')
    plt.ylabel('v')
    plt.title('psf, "%s"'%comment)
    if filename is not None: plt.savefig(filename)

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
    parser.add_option('--time',
                      action='store',
                      dest='ntimeslots',
                      type=int,
                      default=300,
                      help='Number of time slots')
    parser.add_option('--dec',
                      action='store',
                      dest='declination',
                      type=float,
                      default=-90,
                      help='Pointing declination in degrees')
    parser.add_option('--sub',
                      action='store',
                      dest='subarray',
                      type=str,
                      default='mkat',
                      help='Name of subarray as defined in config')
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

    delta_pos=[]
    for ant in numpy.sort(mkat.get_sub(opts.subarray)):
        [delta_N, delta_E, delta_H] = mkat.antennas[ant]
        delta_pos.append([delta_N, delta_E])
    P=numpy.array(delta_pos,dtype='float')
    P = P/l # baseline
    mx = numpy.abs(P).max()

    # number of time slots
    ntimeslots = opts.ntimeslots
    # hour angle range in hours
    ha = np.linspace(-4, 4, ntimeslots)*np.pi/12;

    # source, obs postions
    #latitute
    latitute = (-30. + 42./60. + 48./3600.)*np.pi/180.

    # declination convert in radian
    dec = opts.declination*np.pi/180.

    # matrix to store baseline lengths and azimuth angles
    B = baseline_angles(P)
    #number of antennas
    na = np.size(P, 0)
    #number pair or baseline
    nbl = na*(na-1)/2

##Show UV tracks
    # Plot the uv-Coverage
    plt.figure(figsize=(20,13), facecolor='white')
    plt.hold(True)
    for i in range(nbl):
        uv = track_uv(ha, B[i, 0], 0., B[i, 1], latitute, dec, ntimeslots)
        plt.plot(uv[:,0], uv[:,1], 'b.', markersize=1)
        plt.plot(-uv[:,0], -uv[:,1], 'r.', markersize=1)
    plt.xlabel('u')
    plt.ylabel('v')
    plt.title('Declination %.2f deg: UV coverage'%opts.declination)
    uv = track_uv(ha, B[nbl-1, 0], 0., B[nbl-1, 1], latitute, dec, ntimeslots)
    mb = 5*np.sqrt((uv**2).sum(1)).max();
    plt.xlim(-mb,mb);
    plt.ylim(-mb,mb);
    plt.axis('equal')
    plt.gca().invert_xaxis()
    plt.hold(False)
    if opts.savegraph and opts.verbose: plt.savefig('%s_uv_coverage_%d.png'%(opts.subarray,int(opts.declination)))

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
    filename = None
    if opts.savegraph: filename = '%s_natural_weighting_%d.png'%(opts.subarray,int(opts.declination))
    plot_mask(mask.T, psf, filename, comment="natural weighting")

##Show UV mask and synthesize beam, natural weighting with Gaussian taper
    gauss_kernel = Gauss (npix, npix/20.)
    psf = np.fft.ifftshift(np.fft.ifft2(gauss_kernel*(mask.T))).real
    # plot the mask
    filename = None
    if opts.savegraph: filename = '%s_tapered_gaussian_weighting_%d.png'%(opts.subarray,int(opts.declination))
    plot_mask(gauss_kernel*mask.T, psf, filename, comment="tapered Gaussian weighting")


    if opts.verbose: plt.show()





# -fin-

