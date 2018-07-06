import sys

sys.path.append('/home/kaurov/bin/python/')



au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]


import h5py
import numpy as np
import radmc3dPy.natconst as nc
import matplotlib.pyplot as plt

from sphtool.sphtool import SPHTool
from sphtool.sphtool import RegularGrid


def read_gizmo_dust_sph(fname=None):
    """
	(modified from sphtool example)
    Function to read the example SPH simulation

    Parameters
    -------
        fname   - str
                  Name of the SPH snapshot file (hdf5) to be read
    Returns
    -------
        A dictionary with the following keys
            * npart - Number of particles
            * x     - Cartesian x coordinate of the particle
            * y     - Cartesian y coordinate of the particle
            * z     - Cartesian z coordinate of the particle
            * pmass - Particle mass
            * h     - Smoothing length
            * vx    - Cartesian x velocity compoent of the particle
            * vy    - Cartesian y velocity compoent of the particle
            * vz    - Cartesian z velocity compoent of the particle

        All of the dictionary keys contain a numpy ndarray with a length of npart
    """
    f = h5py.File(fname, 'r')
    return {
        'npart': f['PartType3']['Coordinates'].shape[0],
        'x': f['PartType3']['Coordinates'][:, 0],
        'y': f['PartType3']['Coordinates'][:, 1],
        'z': f['PartType3']['Coordinates'][:, 2],
        'pmass': f['PartType3']['Masses'][:],
        'h': f['PartType3']['Coordinates'][:, 0] * 00 + 0.005,
        'rho': f['PartType3']['Masses'][:] / (4. / 3. * np.pi * f['PartType3']['GrainSize'][:] ** 3),
        'vx': f['PartType3']['Velocities'][:, 0],
        'vy': f['PartType3']['Velocities'][:, 1],
        'vz': f['PartType3']['Velocities'][:, 2],
        'grainsize': f['PartType3']['GrainSize'][:]}


fname = 'snapshot_100.hdf5'
sph = read_gizmo_dust_sph(fname=fname)

temp = sph['z'].copy()
temp[np.where((temp>0.25) & (temp<0.5))[0][::2]] -= 0.25
temp[np.where((temp<0.75) & (temp>0.5))[0][::2]] += 0.25
sph['z'] = temp

# Setting the largest size to be 1 micron
sph['grainsize'] /= (sph['grainsize']/2.).max()
sph['grainsize'] *= 1e-4 * 3e8

sph['rho'] = sph['rho']*0.0 + 2.0
sph['pmass'] = 4./3.*np.pi*sph['rho']*sph['grainsize']**3

# We need to add vector quantities as a three
# dimensional numpy.ndarray objects with the
# first index is the particle indes, second is
# the variable index (here we have only one variable, the velocity),
# and the third is the spatial coordinate index:
vel = np.zeros([sph['x'].shape[0], 1, 3], dtype=np.float64)
vel[:, 0, 0] = sph['vx']
vel[:, 0, 1] = sph['vy']
vel[:, 0, 2] = sph['vz']

bsize = rs*100.

tool = SPHTool(maxSPHParticlesPerCell=10000, maxSPHTreeDepth=10, nThreads=16)
tool.setSPHData(x=((sph['z'] - 0.5)*bsize).astype(np.float64),
                y=((sph['x'] - 0.5)*bsize).astype(np.float64),
                z=((sph['y'] - 0.5)*bsize).astype(np.float64),
                h=(sph['h']*bsize).astype(np.float64),
                rho=sph['rho'].astype(np.float64),
                pmass=sph['pmass'].astype(np.float64),
                vectors=vel)
tool.init()


crd_sys = 'car'
nx = 100
ny = 100
nz = 100
xbound = [-bsize/2., bsize/2.]
ybound = [-bsize/2., bsize/2.]
zbound = [-bsize/2., bsize/2.]

grid = RegularGrid(crd_sys=crd_sys,
                   nx=nx,
                   ny=ny,
                   nz=nz,
                   xbound=xbound,
                   ybound=ybound,
                   zbound=zbound)

scalarFloor = np.zeros(1, dtype=float)+1e-30

tool.regridToRegularGrid(rgrid=grid, scalarFloor=scalarFloor)
tool.writeGrid()
tool.writeScalar(fname="dust_density.binp", ivar=0, nscal=1, binary=True)


### Examine data
fig = plt.figure()
tool.plotScalarSlice2D(incl=0., phi=0.,
                       iscalar=0, npix=400,
                       imsize=bsize, scale='log',
                    cblabel=r'g/cm$^3$')
plt.show()


# fig = plt.figure()
# tool.plotScalarSlice2D(incl=0., phi=0.,
#                        iscalar=0, npix=400,
#                        imsize=1*bsize, cblabel=r'g/cm$^3$')
# plt.show()

### Project on a regular grid:

tool.regridToAMRGrid(maxAMRParticlesPerCell=10000, maxAMRTreeDepth=10, scalarFloor=scalarFloor)
tool.writeGrid()
tool.writeScalar(fname="dust_density.binp", ivar=0, nscal=1, binary=True)










xi       = np.linspace(xbound[0],xbound[1],nx+1)
yi       = np.linspace(ybound[0],ybound[1],ny+1)
zi       = np.linspace(zbound[0],zbound[1],nz+1)
xc       = 0.5 * ( xi[0:nx] + xi[1:nx+1] )
yc       = 0.5 * ( yi[0:ny] + yi[1:ny+1] )
zc       = 0.5 * ( zi[0:nz] + zi[1:nz+1] )
#
# Make the grid
#
qq       = np.meshgrid(xc,yc,zc,indexing='ij')
xx       = qq[0]
yy       = qq[1]
zz       = qq[2]
#

H, edges = np.histogramdd([(sph['z'] - 0.5)*bsize, (sph['x'] - 0.5)*bsize, (sph['y'] - 0.5)*bsize], bins = [xi, yi, zi])
rhod = H*1e-17

with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('0\n')                       # Coordinate system
    f.write('0\n')                       # gridinfo
    f.write('1 1 1\n')                   # Include x,y,z coordinate
    f.write('%d %d %d\n'%(nx,ny,nz))     # Size of grid
    np.savetxt(f,xi.T,fmt=['%13.6e'])    # X coordinates (cell walls)
    np.savetxt(f,yi.T,fmt=['%13.6e'])    # Y coordinates (cell walls)
    np.savetxt(f,zi.T,fmt=['%13.6e'])    # Z coordinates (cell walls)

#
# Write the density file
#

with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    f.write('1\n')                       # Nr of dust species
    data = rhod.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    np.savetxt(f,data.T,fmt=['%13.6e'])  # The data
