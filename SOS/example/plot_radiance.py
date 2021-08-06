# \author Meng Gao
# \date 05/08/2014

import matplotlib.pylab as plt
import numpy as np

#import data
dat1=np.genfromtxt("./radiance.dat")
print dat1.shape
plt.rcParams.update({'font.size': 16, 'font.family': 'serif'})

# create figure and axes
fig=plt.figure()
#plt.ylim([-2,2])
ax=fig.add_subplot(111)

x=dat1[1:-1,0]
y=dat1[1:-1,1]

ax.plot(x,y,'ro')
ax.set_xlabel('Scattering angle ($^\circ$)')
ax.set_ylabel('Diffuse radiance')

# save
plt.savefig('radiance.png', format='png', dpi=200)
