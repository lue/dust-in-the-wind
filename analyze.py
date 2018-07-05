import sys
sys.path.append('/home/kaurov/bin/python/')

from radmc3dPy import *
import matplotlib.pyplot as plt

im = image.readImage()
image.plotImage(im, au=True, log=True, cmap=plt.cm.gist_heat)

analyze.writeDefaultParfile('ppdisk')
par = analyze.readParams()
par.printPar()