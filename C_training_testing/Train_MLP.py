import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.integrate import odeint

import os
cwd=os.getcwd()
from scipy.io import loadmat

x = loadmat(cwd+'/operater_observe.mat')
basis16=x.get('basis4')

x = loadmat(cwd+'/T1.mat')
T1=x.get('T1')

x = loadmat(cwd+'/T2.mat')
T2=x.get('T2')

for s in range(1,2):

    x = loadmat(cwd+'/Coe_'+str(s)+'.mat')
    c_time=x.get('c_time')
    H_Matrix=x.get('H_Matrix')
    
    def S_evolution(y,t):
    
        c=np.zeros([16,1],'complex')
        for i in range(0,16):
            c[i,0]=y[i*2]+1j*y[i*2+1]
    
        dcdt=np.zeros([16,1],'complex')
        dcdt=dcdt+1/1j*np.dot(H_Matrix[0,:,:],c)/5
        dcdt=dcdt+1/1j*np.sin(c_time[0,1]*t+2*np.pi*c_time[0,2])*np.dot(H_Matrix[1,:,:],c)/5
    
        dydt = np.zeros(32)
        for i in range(0,16):
            dydt[i*2]=np.real(dcdt[i,0])
            dydt[i*2+1]=np.imag(dcdt[i,0])
        return dydt

    N=300
    M=100
    t1 = np.linspace(0.01, 5, M*5)
    dt1=t1[1]-t1[0]
    t = np.linspace(0.05, 5, M)
    dt=t[1]-t[0]
    
    x_train=np.zeros([N*M,3])
    y_train=np.zeros([N*M,3])
    t_train=np.zeros([N*M,1])
    IC=np.zeros([N,32])
    for i in range(0,N):
        
        x1_train=np.zeros([M*5,3])
        
        r0=np.zeros(16)
        theta0=np.zeros(16)
        for k in range(0,16):
            r0[k]=random.random()
            theta0[k]=2*np.pi*random.random()
        
        r0=r0/np.sum(r0)
    
        o=0
        for k in range(0,16):
            IC[i,o]=np.sqrt(r0[k])*np.cos(theta0[k])
            o=o+1
            IC[i,o]=np.sqrt(r0[k])*np.sin(theta0[k])
            o=o+1
        
        sol = odeint(S_evolution, IC[i,:], t1)
       
        for k1 in range(0,M*5):        
        
            psi=np.zeros([16,1],'complex')
            for k2 in range(0,16):
                psi[k2]=sol[k1,2*k2]+1j*sol[k1,2*k2+1]
            
            for k2 in range(0,3):
                x1_train[k1,k2]=np.real(np.sum(np.conjugate(psi)*np.dot(basis16[k2,:,:],psi)))
        

        for k1 in range(0,3):
            for k2 in range(0,M*5):
                
                k3=int(k2//5)
                
                if k2%5 == 0:
                    x_train[i*M+k3,k1]=x1_train[k2,k1]
                    
            y_train[i*M+1:(i+1)*M,k1]=(x_train[i*M+1:(i+1)*M,k1]-x_train[i*M:(i+1)*M-1,k1])/dt
            y_train[i*M,k1]=y_train[i*M+1,k1]
        
        t_train[i*M:(i+1)*M,0]=t
        print('ODE',str(s),i)
    
    ######################################    
    # Construct matrix equation Tx=b
    ######################################

    Operator_gather=np.zeros([3,M,256])

    for i in range(0,M):
    
        T=np.zeros([N,256])
        b_gather=np.zeros([N,3])
        for j in range(0,N):
            for k in range(0,256):
                T[j,k]=np.dot(IC[j,:],np.dot(T1[k,:,:],np.transpose(IC[j,:])))
                for k in range(0,3):
                    b_gather[j,k]=x_train[j*M+i,k]
    
        for j in range(0,3):    
            x = np.dot(np.linalg.pinv(T),b_gather[:,j])
            Operator_gather[j,i,:]=x[:]
        print('T1',str(s),i)
    
    ######################################    
    # Construct weight
    ######################################
        
    Weight_G=np.zeros([N*M,3,256])

    for i in range(0,N):
    
        Psi=np.transpose(IC[i,:])
        Psi2=np.dot(np.dot(T2,Psi),Psi)
        for j in range(0,M):
            for k in range(0,3):
        
                o=np.transpose(Operator_gather[k,j,:])
                Weight_G[i*M+j,k,:]=np.dot(Psi2,o)
        print('T2',str(s),i)
        
        ######################################
        # Build model
        ######################################

    from keras.layers import Dense, Input
    from keras import Model
    import keras.backend as K
    from functools import partial
    import tensorflow as tf

    def sample_loss(y_true, y_pred, Weight_G):
    
        Weight=Weight_G[:,0,:]
        h_pred=tf.keras.backend.sum(Weight*y_pred,axis=1)
        Weight=Weight_G[:,1,:]
        x=tf.keras.backend.sum(Weight*y_pred,axis=1)
        h_pred=tf.stack([h_pred, x], axis=1)
        for i in range(2,3):
            Weight1=Weight_G[:,i,:]
            x=tf.keras.backend.sum(Weight1*y_pred,axis=1)
            x=tf.expand_dims(x,axis=1)
            h_pred=tf.concat([h_pred, x], axis=1)

        return K.mean(K.abs(h_pred-y_true)**2,axis=1)

    x = Input(shape=(1,))
    y_true = Input(shape=(3,))
    Weight_G_Input= Input(shape=(3,256,))

    y_pred = Dense(200, activation='tanh')(x)
    y_pred = Dense(200, activation='tanh')(y_pred)
    y_pred = Dense(200, activation='tanh')(y_pred)
    y_pred = Dense(256, activation='tanh')(y_pred)
    model = Model( inputs=[x, y_true, Weight_G_Input], outputs=y_pred)
    model.add_loss( sample_loss( y_true, y_pred, Weight_G_Input) )
    model.compile( loss=None, optimizer='adam' )

    model.fit([t_train, y_train, Weight_G], epochs=100, batch_size=128)

    test_model = Model( inputs=x, outputs=y_pred)

    ######################################
    ## Test
    ######################################
    import scipy.io

    M0=10000
    t_test=np.linspace(0, 5, M0)
    h_pred = test_model.predict(t_test)

    scipy.io.savemat(cwd+'/t_test.mat', mdict={'t_test': t_test})
    scipy.io.savemat(cwd+'/h_pred_'+str(s)+'.mat', mdict={'h_pred': h_pred})