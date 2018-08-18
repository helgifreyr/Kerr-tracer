import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from random import randint

from numpy import loadtxt, cos, sin, linspace, pi, outer, ones, size
from PIL import Image, ImageChops, ImageOps
from papirus import PapirusImage

def draw3D(M,a,r1):
    s, xt, t, xr, r, xth, th, xph, ph = loadtxt('output.dat', unpack=True)
    fline=open('output.dat').readline().rstrip().split()
    rs = float(fline[1])
    M = float(fline[2])
    a = float(fline[3])
    r1 = float(fline[4])
    th1 = float(th[1])
    ph1 = float(ph[1])
    rmax = max(r)

    fig = plt.figure(figsize=(12.8*0.5, 8*0.5))
    ax = fig.gca(projection='3d')
    ax.grid(False)
    ax.axis('off')
    ax.autoscale(enable=False,axis='both')  #you will need this line to change the Z-axis
    ax.set_xbound(-0.3*r1*1.6, 0.3*r1*1.6)
    ax.set_ybound(-0.3*r1*1.6, 0.3*r1*1.6)
    ax.set_zbound(-0.3*r1, 0.3*r1)
    x = r * cos(th) * sin(ph)
    y = r * sin(th) * sin(ph)
    z = r * cos(ph)
    
    ax.plot(x, y, z, color='black', linewidth=1.5)
    
    X,Y,Z = generateSphere(rs,[0,0,0])
    ax.plot_surface(X, Y, Z,  color='black', linewidth=0, alpha=1)
    ax.view_init(elev=randint(-90,90), azim=randint(0,360))
    ax.view_init(elev=0, azim=0)
    plt.savefig('output.png',dpi=800)

def invertImage():
    image = Image.open('output-trimmed.png')
    image = image.convert('L')
    inverted_image = ImageOps.invert(image)
    inverted_image = inverted_image.convert('1')
    inverted_image.save('output-inverted.png')


def trim(file):
    img = Image.open(file+'.png')
    w, h = img.size
    r = 1280./800
    x = int((h-w/r)/2)
    final_img = img.crop((0, x, w, h-x))
    final_img.save(file+'-trimmed.png')

def drawOnScreen():
    image = PapirusImage()
    image.write('output-inverted.png')

def generateSphere(r, c):
    u = linspace(0, 2 * pi, 200)
    v = linspace(0, pi, 200)
    x_c, y_c, z_c = c
    X = r * outer(cos(u), sin(v)) + x_c
    Y = r * outer(sin(u), sin(v)) + y_c 
    Z = r * outer(ones(size(u)), cos(v)) + z_c
    return X,Y,Z

M=1
a=0.5
r1=10
# drawBitmap(M,a,r1)
draw3D(M,a,r1)
trim('output')
invertImage()
# drawOnScreen()
