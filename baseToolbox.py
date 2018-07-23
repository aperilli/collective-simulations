import numpy as np 
import math




def distPointToPoint(x,theta,periodic,L):
    distance = np.zeros((len(x[:,0]),len(x[:,0])))
    X = []
    for k in range(0,len(x[0,:])):
        X1 = np.tile(x[:,k],(len(x[:,k]),1))
        X2 = np.tile(x[:,k],(len(x[:,k]),1)).T
        if k==0 and periodic==1:
            dis1=np.minimum(np.abs(X2-X1),np.abs(L-np.abs(X2-X1)))
            distance += dis1**2
            M1=np.sign(X2-X1)
            O1=np.ones((np.size(M1,0),np.size(M1,1)))
            O1[np.where(dis1!=np.abs(X2-X1))]=-1
            M1=np.multiply(M1,O1)
            dis2=np.minimum(np.abs(X2-X1),np.abs(L-np.abs(X2-X1)))*M1
            X.append(dis2)

        else:
            distance += (X2-X1)**2
            X.append((X2-X1))    
                
    X1 = X[0]*np.cos(-theta)+X[2]*np.sin(-theta)
    X2 = -X[0]*np.sin(-theta)+X[2]*np.cos(-theta)
    angle = np.arctan2(X1,X2)
    distance = np.sqrt(distance)

    np.fill_diagonal(distance,0)
    np.fill_diagonal(angle,0)

    return distance, angle



def angleDifference(A1, A2):
    A = A1 - A2
    A = (A + math.pi) % (2 * math.pi) - math.pi
    return A



def writeCsvRoots(X,name,path,writeMode = 0):
    if writeMode == 0:
        fd = open(path+name,'wb')
    elif writeMode == 1:
        fd = open(path+name,'ab')
    
    wri = X



    np.savetxt(fd,np.c_[wri],delimiter=',',fmt= '%5.5f')
    fd.close()
