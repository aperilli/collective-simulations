#%%%
import matplotlib.animation as animation
from matplotlib import colors as mcolors
import plantsKin_periodic_skeletsize as pk
from math import pi
import numpy as np
import matplotlib.pyplot as plt
import imp
import time
import baseToolbox as bt


LEN=20          #how many time steps in total
Num_roots=21     #should be odd, for symmetry
sig_ang=10.0      #Preceptive noise angular width
sig_P_noise=sig_ang*pi/180
T_P_noise = 10    #Preceptive noise changes every N timesteps  
T_G_noise = 10   #Growth rate noise changes every N timesteps
sig_G_noise = 1.0 #Growth rate noise width
ls=1.0          #separation length between roots (always=1)

color1=(0.7,1,1)
ds_max=0.01     #the distance above which a skeleton will break in two
roots = pk.Roots(Num_roots,ls,T_P_noise,sig_P_noise,T_G_noise,sig_G_noise,theta0 =pi,ds_max=0.01,growth='Apical')

lr=0.3            #repulsion length
la=2*lr         #attraction length
Cr=1.5          #repulsion intensity
Ca=1            #attraction intensity
Var= Cr/Ca * (lr/la)**2 #a figure of merit wrt the H-stable criterion of the Morse potential
print("Cl^2 =" + str(Var))
roots.addCollectiveInteraction(name ='Apical_Morse',repulsionZone=lr,attractionZone=la,repulsionIntensity=Cr,attractionIntensity=Ca)
#adding other interactions
#roots.addCollectiveInteraction(name ='Apical')
#roots.addInteractions(name = 'Tropism' ,intensity=-0.5,direction = -pi)
#roots.addInteractions(name = 'Proprioception' ,intensity=-0.25)

separation=10    #marker separation steps
inrun_timestep_separation = 10 #add a marker every N timesteps
Markers_length=int(1/ds_max/separation-1)
for root in roots.roots:
    root.set_markers(separation)

roots.update()  #first timestep        

#animation stuff
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)


#plotting
fig = plt.figure(figsize=(10.80, 10.80), dpi=100)
fig.tight_layout()
ax = fig.add_subplot(111)
ax.axis('off')
#plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
lines = []
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
#resent the lines for the animation 
for indind in range(0,5):   #more than one copy of lines, for periodic line crossings
    for root in roots.roots:
        obj1=ax.plot([],[], antialiased=True)[0]
        plt.setp(obj1, color=mcolors.hsv_to_rgb(color1), linewidth=2)
        lines.append(obj1)
for root in roots.roots:    #extra lines for the markers
        obj1=ax.plot([],[], 'b.')[0]
        lines.append(obj1)

    
ax.set_xlim(0,Num_roots)
ax.set_ylim(-20,0.2)
ax.set_aspect('equal', 'box')
    
def init():         #lines initialization for the animation
    for line in lines:
        line.set_data([],[])
    return lines

def animate(i):     #an animation timestep
    t0=time.time()
    roots.update()
    t1=time.time()
    ax.set_xlim(0,Num_roots)
    ax.set_ylim(-20,0.2)
    ax.set_aspect('equal', 'box')
#    print(-1/(t0-t1))
    t0=t1
    ind_lines=0
    for k in range(0,len(roots.roots)):
        DIFF=np.diff(roots.roots[k].x[:,0])
        if np.size(DIFF[np.where(np.abs(DIFF)>(Num_roots-3))]):     #identifing PBC crosssings
            ind1=np.argwhere((np.abs(DIFF)>(Num_roots-3)))
            nni=1
            nnf=-1
            ind1=np.insert(ind1,0,-nni)
            for jnd1 in range(np.size(ind1)-1):                     #splitting the non-continuous lines
                ind_lines+=1
                lines[len(roots.roots)+ind_lines].set_data(roots.roots[k].x[ind1[jnd1]+nni:ind1[jnd1+1]+nnf,0],roots.roots[k].x[ind1[jnd1]+nni:ind1[jnd1+1]+nnf,2])
            lines[k].set_data(roots.roots[k].x[ind1[-1]+nni::,0],roots.roots[k].x[ind1[-1]+nni::,2])
        else:
            lines[k].set_data(roots.roots[k].x[:,0],roots.roots[k].x[:,2])
        #markers:    
        Markers=np.zeros((1,3))
        q=0
        if np.mod(i,inrun_timestep_separation)==0:
            for root in roots.roots:
                root.skeleton[-2].set_marker()
        for skel in roots.roots[k].skeleton:
            if skel.marker==1:
                if q==0:
                    for ind22 in range(0,3):
                        Markers[0,ind22]=np.copy(skel.x[ind22])
#                    print(Markers)
#                    print(np.shape(Markers))
                else:
                    dumm=np.zeros((1,3))
                    for ind22 in range(0,3):
                        dumm[0,ind22]=np.copy(skel.x[ind22])
                    Markers=np.append(Markers,np.copy(dumm),axis=0)
                q+=1
#        print(Markers)
        lines[-k-1].set_data(Markers[:,0],Markers[:,2])
    return lines


#saving animation
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=LEN, interval=30, blit=True)
anim.save('basic_animation.mp4', fps=10, extra_args=['-vcodec', 'libx264'])
#%%
fig = plt.figure(figsize=(10.80, 10.80), dpi=100)
fig.tight_layout()
ax = fig.add_subplot(111)
ax.axis('off')

#k=10
for k in range(0,5):
    X=lines[k].get_xdata()
    Y=lines[k].get_ydata()
    ax.plot(X,Y)
#for line in lines:
#    ax.plot(line.get_xdata(),line.get_ydata())
#    X=line.get_xdata()
#    print(X[0:5])
