# \author Meng Gao
# \date April 30, 2014
# \brief This code is used to plot the photo path from the ray tracing program

import matplotlib.pyplot as plt
import numpy as np

# create figure and axes
fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot(111, aspect='equal')

# axis
ax.axis('off')
plt.axis([-2, 1.5, -1.5, 1.5])

# plot a disk
circle2=plt.Circle((.0,.0),1.0,color='b', alpha=0.1)
ax.add_artist(circle2)

#load data for line plot
dat2=np.genfromtxt("./path.dat")
d2=len(dat2)
# plot inner rays
xv=dat2[:,1]
yv=dat2[:,2]
rays = plt.Line2D(
  xv,yv, 
    linewidth=1, 
    color='r'
    )
ax.add_artist(rays)

# plot incident rays
xvi=[-2,xv[0]]
yvi=[yv[0],yv[0]]
rayi = plt.Line2D(
  xvi,yvi, 
    linewidth=3, 
    color='r'
    )
ax.add_artist(rayi)

for i in range(0,d2):
    # transmitted rays
    ax.arrow(dat2[i,1], dat2[i,2], 0.15*dat2[i,3],0.15*dat2[i,4], # x0, y0, dx, dy
            head_width=0.05,
            alpha=0.9,
            color='r'
            )
    # reflected rays
    ax.arrow(dat2[i,1], dat2[i,2], 0.15*dat2[i,5], 0.15*dat2[i,6], 
            head_width=0.05,
            alpha=1,
            color='r'
            #head_length=0.1
            #, fc='k', ec='k'
            )

# save
plt.savefig('ray_path.png', format='png', dpi=200)
