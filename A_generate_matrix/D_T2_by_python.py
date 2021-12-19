import numpy as np

N=2**4
X1=np.zeros([1,N,2*N],'complex')
X3=np.zeros([N,2*N,1],'complex')
for i in range(0,N):
    X1[0,i,2*i]=1
    X1[0,i,2*i]=-1j
    X3[i,2*i-1,0]=1
    X3[i,2*i,0]=1j

X2=np.zeros([N,N,N**2],'complex')

o=0
for i in range(0,N):
    for j in range(i,N):
        
        if i==j :
            X2[i,j,o]=1
            o=o+1
        else:
            X2[i,j,o]=1
            X2[i,j,o+1]=1j
            o=o+2


for i in range(0,N):
    for j in range(0,i-1):
        for k in range(0,N**2):
            X2[i,j,k]=np.conj(X2[j,i,k]);

T2=np.zeros([N**2,N**2,2*N,2*N],'complex')
for i1 in range(0,N**2):
    for i2 in range(0,N**2):
        for i3 in range(0,2*N):
            for i4 in range(0,2*N):
            
                x1=X1[:,:,i3]
                x1=np.reshape(x1,(1,N))
                x2=X2[:,:,i1];
                x2=np.reshape(x2,(N,N))
                x3=X2[:,:,i2]
                x3=np.reshape(x3,(N,N))
                x4=X3[:,i4,:]
                x4=np.reshape(x4,(N,1))
            
                T2[i1,i2,i3,i4]=1j*np.dot(np.dot(x1,(np.dot(x2,x3)-np.dot(x3,x2))),x4)
            

    print('1',i1/N**2)


for i1 in range(0,N**2):
    for i1 in range(0,N**2):
        for i3 in range(0,2*N):
            for i4 in range(0,i3):
                if (T2[i1,i2,i3,i4]==-T2[i1,i2,i4,i3] and T2[i1,i2,i3,i4]!=0):
                    T2[i1,i2,i3,i4]=0
                    T2[i1,i2,i4,i3]=0

    print('2',i1/N**2)


#save('T2.mat','T2')