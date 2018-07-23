import numpy as np
import math
import time
import baseToolbox as bt
from math import pi
import copy

InteractionsName = ['Tropism','ApicalTropism','Proprioception','ApicalRelative']
CollectiveInteractionsName = ['Apical','Global','Apical_Morse']
GrowthMode = ['Apical','Exponential']



class skeletonElements:
    def __init__(self,x,theta0,curvature0,s,ds,psi0 = 0,psiC0 = 0,psiG0 = 0,dt = .1,L=3,growth = {'name':'no','intensity':1,'direction':0},T_G_noise=1,sig_G_noise=0):
        
        self.x = x
        self.theta = theta0
        self.curvature = curvature0
        self.ds = ds
        self.s = s
        self.psiC = psiC0
        self.psiG = psiG0
        self.psi0 = psi0
        self.dt = dt
        self.growth = growth.copy()
        self.interactions = []
        self.thetaApical = theta0
        self.boxsize=L
        self.istoolong=0
        self.marker = 0
        self.sig_G_noise=sig_G_noise
        self.T_G_noise=T_G_noise
        self.G_noise=float(0)
        
    def updateCurvilinearAbscissa(self,s):
        self.s = s

    def updateApicalAngle(self,theta):
        self.thetaApical = theta

    def updateOrientation(self,theta):
        self.theta = theta  

    def updateSpatialPosition(self,x0,theta0,psi0,ds0):
        self.x[0]  = np.mod(x0[0] + ds0 * (np.sin(theta0) * np.cos(psi0)),self.boxsize)
        self.x[1]  = x0[1] + ds0 * (np.sin(theta0) * np.sin(psi0))
        self.x[2]  = x0[2] + ds0 * (np.cos(theta0))

    def istoolong_update(self,ds_max):
        if self.ds>=2*ds_max:
            self.istoolong=1
        else:
            self.istoolong=0
            
    def set_marker(self):
        self.marker = 1
        
    def addInteractions(self,name,intensity=0,direction = 0):
        
        self.interactions.append({'name':name,'intensity':intensity,'direction':direction})
        
    def update(self,time_counter):
        t0 = time.time()

        delta = 0
        deltaParrallel = 0
        deltaPerpendicular = 0
        for interaction in self.interactions:
            #print('interaction : '+str(interaction['direction']))
            interactionName = interaction['name']
            if interaction['name'] == 'Proprioception':
                deltaParrallel +=interaction['intensity'] * self.curvature 
            else:
                if interaction['name'] == 'Tropism':
                    delta += interaction['intensity'] * np.sin(bt.angleDifference(self.theta,interaction['direction'])) #sin is added!!!!!!
                elif interaction['name'] == 'ApicalTropism':
                    delta += interaction['intensity'] * np.sin(bt.angleDifference(self.thetaApical,interaction['direction'])) #sin is added!!!!!!
                elif interaction['name'] == 'ApicalRelative':
                    delta += interaction['intensity'] * np.sin(interaction['direction']) #sin is added!!!!!!
        #print('angle : '+str(self.thetaApical))
        ###################################################         
        #DELTA MUST BE SMALLER THAN 1. MUST FIX INTESITIES#
        ####################################################
        deltaParrallel +=   delta * np.cos(self.psiG - self.psiC)
        deltaPerpendicular +=   delta * np.sin(self.psiG - self.psiC)
        
        
        self.curvature += deltaParrallel * self.growth['growthRate'] * self.dt 
        if np.abs(self.curvature) >0:
            self.psiC = (deltaPerpendicular * self.growth['growthRate'] * self.dt)/self.curvature
        if self.growth['name'] in GrowthMode:
            if self.growth['growthRate']:
                if np.mod(time_counter,self.T_G_noise)==1:
                    self.G_noise= float(np.random.uniform(-self.sig_G_noise,self.sig_G_noise,1))
                self.ds += self.ds * (self.growth['growthRate'] +self.G_noise)* self.dt 


class Plant:
    def __init__(self,x0,theta0=0,curvature0=0,length0 = 1,psi0 = 0,psiC0 = 0,psiG0 = 0,ds_max=0.01,dt=.1,L=15,growth = 'no',growthZone = 1,growthRate = 1,T_G_noise=1,sig_G_noise=0):
        
        self.x = []
        self.s = []
        self.theta = []
        self.curvature = [] 
        
        self.boxsize=L
        self.x0 = x0
        self.theta0 = theta0
        self.curvature0 = curvature0
        
        self.length0 = length0
        self.length = length0
        self.ds = ds_max
        self.ds_max=ds_max
        self.psi0 = psi0
        self.psiG0 = psiG0
        self.psiC0 = psiC0
        self.dt = dt
        self.skeleton = []
        self.skeleton_copy=[]
        
        self.interactions = []
        x = self.x0
        theta = self.theta0
        curvature = self.curvature0

        self.growth = {'name':growth,'growthZone':growthZone,'growthRate':growthRate}
        self.T_G_noise=T_G_noise
        self.sig_G_noise=sig_G_noise
        
        s = 0

        self.skeleton.append(skeletonElements(np.copy(x),theta,curvature,0,self.ds,self.psi0,psiC0,psiG0,self.dt,self.boxsize,growth=self.growth,T_G_noise=self.T_G_noise,sig_G_noise=self.sig_G_noise))
        
        
        for k in range(1,int(np.floor(self.length/self.ds))):
            
            x[0]  = np.mod(x0[0] + self.ds * (np.sin(theta0) * np.cos(psi0)),self.boxsize)
            x[1]  = x[1] + self.ds * (np.sin(theta0) * np.sin(psi0))
            x[2]  = x[2] + self.ds * (np.cos(theta0))
            s += self.ds
            self.skeleton.append(skeletonElements(np.copy(x),theta,curvature,s,self.ds,self.psi0,psiC0,psiG0,self.dt,self.boxsize,growth=self.growth,T_G_noise=self.T_G_noise,sig_G_noise=self.sig_G_noise))
        self.flatten()

    def set_markers(self,separation):
        for k in range(1,len(self.skeleton)):
            if np.mod(k,separation)==0:
                self.skeleton[k].set_marker()
                
        
    def updateSpatialPosition(self):
        s = 0
        for k in range(1,len(self.skeleton)):
            s += self.skeleton[k-1].ds
            self.skeleton[k].updateCurvilinearAbscissa(s)
            self.skeleton[k].updateOrientation(bt.angleDifference(self.skeleton[k-1].theta, - self.skeleton[k-1].curvature * self.skeleton[k-1].ds))
            self.skeleton[k].updateSpatialPosition(self.skeleton[k-1].x,self.skeleton[k-1].theta,self.skeleton[k-1].psiC,self.skeleton[k-1].ds)
        for skel in self.skeleton:
            skel.updateApicalAngle(self.skeleton[-1].theta)


    def addInteractions(self,name,intensity=0,direction = 0):
        if name in InteractionsName:
            self.interactions.append({'name':name,'intensity':intensity,'direction':direction})
            for skel in self.skeleton:
                skel.addInteractions(name,intensity,direction)            
        else:
            print(' --- '+name+' is not part of the known interactions ')
            print(' --- please use one of the following interaction :')
            for names in InteractionsName:
                print(' --- --- '+str(names))

    def flatten(self):
        
        self.x=np.array([skel.x for skel in self.skeleton])
        self.s=np.array([skel.s for skel in self.skeleton])
        self.theta=np.array([skel.theta for skel in self.skeleton])
        self.curvature=np.array([skel.curvature for skel in self.skeleton])
        self.length = self.s[-1]

    def break_skeleton(self):
        self.skeleton_copy=[]        
        for skel in self.skeleton:
            skel.istoolong_update(self.ds_max)
            if skel.istoolong==0:
                self.skeleton_copy.append(skeletonElements(skel.x,skel.theta,skel.curvature,skel.s,skel.ds,skel.psi0,skel.psiC,skel.psiG,skel.dt,skel.boxsize,growth=skel.growth,T_G_noise=skel.T_G_noise,sig_G_noise=skel.sig_G_noise))
                if skel.marker==1:
                    self.skeleton_copy[-1].set_marker()
            else:
                DS_temp=skel.ds/2
                self.skeleton_copy.append(skeletonElements(skel.x,skel.theta,skel.curvature,skel.s,DS_temp,skel.psi0,skel.psiC,skel.psiG,skel.dt,skel.boxsize,growth=skel.growth,T_G_noise=skel.T_G_noise,sig_G_noise=skel.sig_G_noise))
                if skel.marker==1:
                    self.skeleton_copy[-1].set_marker()
                X_temp=np.copy(skel.x)
                X_temp[0]  = np.mod(X_temp[0] + DS_temp * (np.sin(skel.theta) * np.cos(skel.psi0)),skel.boxsize)
                X_temp[1]  = X_temp[1] + DS_temp * (np.sin(skel.theta) * np.sin(skel.psi0))
                X_temp[2]  = X_temp[2] + DS_temp * (np.cos(skel.theta))
                self.skeleton_copy.append(skeletonElements(np.copy(X_temp),skel.theta,skel.curvature,skel.s+DS_temp,DS_temp,skel.psi0,skel.psiC,skel.psiG,skel.dt,skel.boxsize,growth=skel.growth,T_G_noise=skel.T_G_noise,sig_G_noise=skel.sig_G_noise))
#                if skel.marker==1:
#                    self.skeleton_copy[-1].set_marker()
        self.skeleton=copy.deepcopy(self.skeleton_copy)
        self.skeleton_copy=[]
        self.updateSpatialPosition()
        self.flatten()

    def updateGrowth(self):
        
        if self.growth['name'] in GrowthMode:
            for skel in self.skeleton:
                
                if self.growth['name'] == 'Exponential':
                    skel.growth['growthRate'] = self.growth['growthRate'] 
                if self.growth['name'] == 'Apical':
                    
                    if (self.length-skel.s) < self.growth['growthZone']:

                        skel.growth['growthRate'] = self.growth['growthRate']
                    else:
                        
                        skel.growth['growthRate'] = 0

    def update(self,time_counter):
#        self.break_skeleton()
        self.updateGrowth()
        
        for skel in self.skeleton:
            skel.update(time_counter)
        self.updateSpatialPosition()

        self.flatten()

    def updateCollectiveInteraction(self):
        for skel in self.skeleton:
            skel.interactions = self.interactions


class Roots:
    def __init__(self,N,dx,T_P_noise,sig_P_noise,T_G_noise,sig_G_noise,theta0,ds_max=0.01,dt=0.1,growth = 'no',growthRate=1,growthZone = 1):
        self.time_counter=0
        self.N = N
        self.dx = dx
        self.sig_P_noise=sig_P_noise
        self.T_P_noise=T_P_noise
        self.sig_G_noise=sig_G_noise
        self.T_G_noise=T_G_noise
        self.roots = []
        self.interactions = []
        self.collectiveInteractions = []
        self.collectiveInteractionsList = [] 
        self.P_Noise_vec=[]
        self.Noise_P_vec=[]
#        self.G_Noise_vec=[]
        for k in range(0,N):
            self.roots.append(Plant(x0 = [k*dx,0,0],theta0 = theta0,ds_max=ds_max,dt=dt,L=N*dx,growth = growth,growthRate=growthRate,growthZone = growthZone,sig_G_noise=self.sig_G_noise,T_G_noise=self.T_G_noise))
        
    def update(self):
        self.time_counter+=1
        if self.collectiveInteractions:

            self.collectiveComputation()
            for k in range(0,self.N):
#                self.roots[k].break_skeleton()               
                self.collectiveInteractionsList[k]['intensity'] = self.intensityCollective[k]
                self.collectiveInteractionsList[k]['direction'] = self.directionCollective[k]
                self.roots[k].break_skeleton()
                self.roots[k].updateCollectiveInteraction()
                self.roots[k].update(self.time_counter)
        else:
            for root in self.roots:
                root.break_skeleton()
                root.update(self.time_counter)

    def addInteractions(self,name,intensity=0,direction = 0):
        if name in InteractionsName:
            self.interactions.append({'name':name,'intensity':intensity,'direction':direction})
            for root in self.roots:
                root.addInteractions(name,intensity,direction)
        

        else:
            print(' --- '+name+' is not part of the known interactions ')
            print(' --- please use one of the following interaction :')
            for names in InteractionsName:
                print(' --- --- '+str(names))

    def addCollectiveInteraction(self,name,repulsionZone,attractionZone,repulsionIntensity,attractionIntensity):
        if name in CollectiveInteractionsName:
            self.collectiveInteractions.append({'name':name,'repulsionZone':repulsionZone,'attractionZone':attractionZone,'repulsionIntensity':repulsionIntensity,'attractionIntensity':attractionIntensity})
            for root in self.roots:
                root.addInteractions('ApicalRelative',0,0)
                self.collectiveInteractionsList.append(root.interactions[-1])
        else:
            print(' --- '+name+' is not part of the known collective interactions ')
            print(' --- please use one of the following collective interaction :')
            for names in CollectiveInteractionsName:
                print(' --- --- '+str(names))

    def collectiveComputation(self):
        self.flatten()
        for interaction in self.collectiveInteractions:
            if interaction['name'] == 'Apical':
                self.tipDistance()
                interactionTip =np.copy(self.distanceTip)
                interactionTip[(self.distanceTip>0) & (self.distanceTip<interaction['repulsionZone'])]=-1
                interactionTip[(self.distanceTip>interaction['repulsionZone']) & (self.distanceTip<interaction['attractionZone'])]=1
                interactionTip[self.distanceTip>interaction['attractionZone']]=0
                
                self.directionCollective = np.sum(interactionTip*self.alphaTip,0)/(self.N-1)
                
                self.intensityCollective = self.directionCollective*0
                self.intensityCollective[np.sum(np.abs(interactionTip),0)>0] =1.0
            
            if interaction['name'] == 'Apical_Morse':
                self.tipDistance()
                interactionTip =np.copy(self.distanceTip)
                R=np.copy(self.distanceTip)
                interactionTip=(interaction['repulsionIntensity']/interaction['repulsionZone'])*np.exp(-R/interaction['repulsionZone']) - (interaction['attractionIntensity']/interaction['attractionZone'])*np.exp(-R/interaction['attractionZone'])
                np.fill_diagonal(interactionTip,0)
                phasors=np.copy(self.distanceTip)
                phasors=-interactionTip*np.exp(1j*self.alphaTip)
                
                self.directionCollective = np.angle(np.sum(phasors,0))
                self.directionCollective[np.where(np.abs(np.sum(phasors,0))<1e-15)]=pi
                if np.mod(self.time_counter,self.T_P_noise)==1:
                    self.P_Noise_vec =np.random.normal(0,self.sig_P_noise, (np.size(self.directionCollective,0)))
                self.directionCollective=self.directionCollective+self.P_Noise_vec
#                self.directionCollective=self.P_Noise_vec #noise check
                self.directionCollective[np.where(np.abs(self.directionCollective)==pi)]=pi
                self.intensityCollective = self.directionCollective*0
                self.intensityCollective[np.sum(np.abs(interactionTip),0)>0] =1.0



    def flatten(self):
        
        self.xTip=np.array([root.x[-1] for root in self.roots])
        self.thetaTip=np.array([root.theta[-1] for root in self.roots])

    def tipDistance(self):
        self.distanceTip,self.alphaTip = bt.distPointToPoint(self.xTip,self.thetaTip,1,(self.N)*self.dx)
