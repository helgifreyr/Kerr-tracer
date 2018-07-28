import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
from random import randint


from numpy import sqrt, pi, sin, cos



def draw3D(M,a,r1):
    s, xt, t, xr, r, xth, th, xph, ph = loadtxt('output.dat', unpack=True)
    rs = M+sqrt(M**2-a**2)
    rs = 2*M

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.grid(False)
    ax.axis('off')
    ax.autoscale(enable=False,axis='both')  #you will need this line to change the Z-axis
    ax.set_xbound(-1.1*r1, 1.1*r1)
    ax.set_ybound(-1.1*r1, 1.1*r1)
    ax.set_zbound(-1.1*r1, 1.1*r1)
    x = r * cos(th) * sin(ph)
    y = r * sin(th) * sin(ph)
    z = r * cos(ph)
    
    ax.plot(x, y, z, color='black', linewidth=0.5)
    
    u = linspace(0, 2 * pi, 200)
    v = linspace(0, pi, 200)

    X = rs * outer(cos(u), sin(v))
    Y = rs * outer(sin(u), sin(v))
    Z = rs * outer(ones(size(u)), cos(v))
    
    ax.plot_surface(X, Y, Z,  color='black', linewidth=0, alpha=1)
    ax.view_init(elev=20, azim=0)
    ax.view_init(elev=randint(-120,120), azim=randint(0,360))
    savefig('output.png',dpi=800)

def invertImage():
    image = Image.open('output.png')
    image = image.convert('L')
    inverted_image = PIL.ImageOps.invert(image)
    inverted_image = inverted_image.convert('1')
    inverted_image.save('output-inverted.png')

    
from pylab import figure, polar, plot, xlabel, ylabel, grid, hold, legend, title, savefig, show, subplot
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from numpy import loadtxt, sin, meshgrid, sqrt, linspace, pi, outer, ones, size
from PIL import Image
import PIL.ImageOps    
M=1
a=0.5
r1=10
#drawBitmap(M,a,r1)
draw3D(M,a,r1)
#trim('output')
invertImage()

