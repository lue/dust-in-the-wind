import sys
sys.path.append('/home/kaurov/bin/python/')

from radmc3dPy import *
import matplotlib.pyplot as plt
import matplotlib.pylab as plb







data = analyze.readData(ddens=False)


c = plb.contourf(data.grid.x, data.grid.y, np.log10(data.rhodust[:,:,0,0].T), 30)
plb.xlabel('r [AU]')
plb.ylabel(r'$\pi/2-\theta$')
# plb.xscale('log')
plt.show()



data = analyze.readData(dtemp=True)
data.readDustTemp()

c = plb.contourf(data.grid.x/natconst.au, data.grid.y/natconst.au, data.dusttemp[:,:,16,0].T, 30)
plb.xlabel('r [AU]')
plb.ylabel(r'$\pi/2-\theta$')
# plb.xscale('log')
cb = plb.colorbar(c)
cb.set_label('T [K]', rotation=270.)
plt.show()

c = plb.contour(data.grid.x/natconst.au, data.grid.y/natconst.au, data.dusttemp[:,:,0,0].T, 10,  colors='k', linestyles='solid')
plb.clabel(c, inline=1, fontsize=10)
plt.show()









im = image.readImage(); image.plotImage(im, au=True, log=True, cmap=plt.cm.gist_heat)













analyze.writeDefaultParfile('ppdisk')
par = analyze.readParams()
par.printPar()




../radmc-3d/version_0.41/src/radmc3d mctherm setthreads 24
