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
        'h': f['PartType3']['Coordinates'][:, 0] * 00 + 0.01,
        'rho': f['PartType3']['Masses'][:] / (4. / 3. * np.pi * f['PartType3']['GrainSize'][:] ** 3),
        'vx': f['PartType3']['Velocities'][:, 0],
        'vy': f['PartType3']['Velocities'][:, 1],
        'vz': f['PartType3']['Velocities'][:, 2],
        'grainsize': f['PartType3']['GrainSize'][:]}


fname = 'snapshot_100.hdf5'
sph = read_gizmo_dust_sph(fname=fname)

# Setting the largest size to be 1 micron
sph['grainsize'] /= (sph['grainsize']/2.).max()
sph['grainsize'] *= 1e-4 * 1e10

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

tool = SPHTool(maxSPHParticlesPerCell=1000, maxSPHTreeDepth=10, nThreads=4)
tool.setSPHData(x=((sph['z'] - 0.5)*bsize).astype(np.float64),
                y=((sph['x'] - 0.5)*bsize).astype(np.float64),
                z=((sph['y'] - 0.5)*bsize).astype(np.float64),
                h=(sph['h']*bsize).astype(np.float64),
                rho=sph['rho'].astype(np.float64),
                pmass=sph['pmass'].astype(np.float64),
                vectors=vel)
tool.init()


### Examine data
fig = plt.figure()
tool.plotScalarSlice2D(incl=0., phi=0.,
                       iscalar=0, npix=400,
                       imsize=bsize, scale='log',
                    cblabel=r'g/cm$^3$')
plt.show()


fig = plt.figure()
tool.plotScalarSlice2D(incl=0., phi=0.,
                       iscalar=0, npix=400,
                       imsize=1*bsize, cblabel=r'g/cm$^3$')
plt.show()

### Project on a regular grid:



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
tool.writeScalar(fname="dust_density.inp", ivar=0, nscal=1, binary=False)

# tool.regridToAMRGrid(maxAMRParticlesPerCell=10, maxAMRTreeDepth=12, scalarFloor=scalarFloor)
# tool.writeGrid()
# tool.writeScalar(fname="dust_density.inp", ivar=0, nscal=1, binary=False)
