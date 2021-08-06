# \author Meng Gao
# \date April 30, 2014
# \brief This code is used to plot the phase function from the ray tracing program

import matplotlib.pylab as plt
import numpy as np

#import data
dat1=np.genfromtxt("./p11.dat")
plt.rcParams.update({'font.size': 16, 'font.family': 'serif'})

# create figure and axes
fig=plt.figure()
ax=fig.add_subplot(111)

x=dat1[1:-1,0]
y=dat1[1:-1,1]

ax.plot(x,y,'r')
ax.set_yscale('log')
ax.set_xlabel('Scattering angle ($^\circ$)')
ax.set_ylabel('Phase function')

# save
plt.savefig('phase_function.png', format='png', dpi=200)
