#
# Import NumPy for array handling
#
import numpy as np
#
# Import plotting libraries (start Python with ipython --matplotlib)
#
#from mpl_toolkits.mplot3d import axes3d
#from matplotlib import pyplot as plt
#
# Some natural constants
#

au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]

#
# Monte Carlo parameters
#

nphot    = 1000000

#
# Grid parameters
#

nx       = 100
ny       = 100
nz       = 100
sizex    = 10.*rs
sizey    = 10.*rs
sizez    = 10.*rs

#
# Model parameters
#
radius   = 1200*rs
rho0     = 1e-13
#
# Star parameters
#
mstar    = ms
rstar    = rs*10.
tstar    = ts
pstar    = np.array([0.,0.,100.*rs])
#
# Make the coordinates
# Write the wavelength_micron.inp file
#
lam1     = 0.1e0
lam2     = 7.0e0
lam3     = 25.e0
lam4     = 1.0e4
n12      = 20
n23      = 100
n34      = 30
lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
lam      = np.concatenate([lam12,lam23,lam34])
nlam     = lam.size

#
# Write the wavelength file
#

with open('wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(nlam))
    np.savetxt(f,lam.T,fmt=['%13.6e'])

#
#
# Write the stars.inp file
#


mstar    = ms
rstar    = rs*10.
tstar    = ts
pstar    = np.array([0.,0.,100.*rs])

with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('1 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
    np.savetxt(f,lam.T,fmt=['%13.6e'])
    f.write('\n%13.6e\n'%(-tstar))

#
# Write the grid file
#
# with open('amr_grid.inp','w+') as f:
#     f.write('1\n')                       # iformat
#     f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
#     f.write('0\n')                       # Coordinate system
#     f.write('0\n')                       # gridinfo
#     f.write('1 1 1\n')                   # Include x,y,z coordinate
#     f.write('%d %d %d\n'%(nx,ny,nz))     # Size of grid
#     np.savetxt(f,xi.T,fmt=['%13.6e'])    # X coordinates (cell walls)
#     np.savetxt(f,yi.T,fmt=['%13.6e'])    # Y coordinates (cell walls)
#     np.savetxt(f,zi.T,fmt=['%13.6e'])    # Z coordinates (cell walls)
# #
# # Write the density file
# #
# with open('dust_density.inp','w+') as f:
#     f.write('1\n')                       # Format number
#     f.write('%d\n'%(nx*ny*nz))           # Nr of cells
#     f.write('1\n')                       # Nr of dust species
#     data = rhod.ravel(order='F')         # Create a 1-D view, fortran-style indexing
#     np.savetxt(f,data.T,fmt=['%13.6e'])  # The data
#
# Dust opacity control file
#

with open('dustopac.inp','w+') as f:
    f.write('2               Format number of this file\n')
    f.write('1               Nr of dust species\n')
    f.write('============================================================================\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('silicate        Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')

#
# Write the radmc3d.inp control file
#

with open('radmc3d.inp','w') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 1000\n')
    f.write('iranfreqmode = 1\n')
    f.write('istar_sphere = 1')

































