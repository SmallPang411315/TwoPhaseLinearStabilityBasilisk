import math
import matplotlib.pyplot as plt
plt.rc('font',family='Times New Roman') 
import numpy as np
import matplotlib.image as mpimg
from scipy.interpolate import griddata
from matplotlib import rcParams
from scipy import interpolate
import scipy.interpolate as spi

fileName          = 'EigenPlot'
readDir           = r'./'
saveDir2          = r'./'

###Font info
font_axis = {'family': 'Times New Roman',
         'style':'italic',
         'weight': 'normal',
         'size': 40,
         }
font_legend = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 40,
         }     

font_eqs = {
    "font.family":'serif',
    "font.size": 40,
    "mathtext.fontset":'stix',
}
ticksize = 40
rcParams.update(font_eqs)
######################################TECPLOT_TIMESLICE#####################################
f=open(r'./EigenVec_205','r')
x=[]#创建三个列表用于存放点的数据
y=[]
h1=[]
h2=[]
for line in f:
    p=line[:]
    a=p.split(' ')
    x.append(float(a[1]))
    y.append(float(a[2]))
    h1.append(float(a[3]))
    h2.append(float(a[4]))
f.close()

points=[]
for i in range(len(x)):
    point=[]
    point.append(x[i])
    point.append(y[i])   
    points.append(point)
points=np.array(points)

xi=np.linspace(min(x),max(x),500)
yi=np.linspace(min(y),max(y),500)

xi,yi=np.meshgrid(xi,yi)#网格化

fig, ax = plt.subplots(figsize=(32, 4))

#plt.subplot(1,2,1)
zi=griddata(points,h1,(xi,yi),method='cubic')#内插
v = np.linspace(-5e-3, 5e-3, 20, endpoint=True)
Power0 = plt.contourf(xi, yi, zi, v,cmap = plt.cm.coolwarm, extend="both")

#plt.contour(Power0,v,colors = 'black',linewidths = .5)
#theta = np.linspace(0, 2 * np.pi, 200)
#x = 0.5*np.cos(theta)
#y = 0.5*np.sin(theta)
#plt.fill(x, y, 'w')

plt.xticks(fontproperties='Times New Roman', size=ticksize)
plt.yticks(fontproperties='Times New Roman', size=ticksize)
plt.xlabel(r'x', font_axis)
plt.ylabel(r'y', font_axis,rotation=360, labelpad=10)
plt.title(label=r'$\hat{u}\prime,We=5$', fontdict=font_legend)

#plt.subplot(1,2,2)
#zi=griddata(points,h2,(xi,yi),method='cubic')#内插
#v = np.linspace(-8e-3, 8e-3, 20, endpoint=True)
#Power0 = plt.contourf(xi, yi, zi, v,cmap = plt.cm.coolwarm, extend="both")

##plt.contour(Power0,v,colors = 'black',linewidths = .5)
##theta = np.linspace(0, 2 * np.pi, 200)
##x = 0.5*np.cos(theta)
##y = 0.5*np.sin(theta)
##plt.fill(x, y, 'w')

#plt.xticks(fontproperties='Times New Roman', size=ticksize)
#plt.yticks(fontproperties='Times New Roman', size=ticksize)
#plt.xlabel(r'x', font_axis)
#plt.title(label=r'$\hat{v}\prime$', fontdict=font_legend)

plt.savefig(saveDir2 + fileName + '.eps', format='eps', dpi=500, bbox_inches='tight')
