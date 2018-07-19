# -*- coding: utf-8 -*-
"""
Created on Sat Jun 16 11:46:10 2018

@author: Fahimi
"""
import numpy as np
from __future__ import division

# data format should be channels by time by trials
# p is model order
data=data*(10**12)
def GeneralLinearKalmanFilter(data,p, sampling_rate,outputname):
    m=np.shape(data)[0] # number of channels
    N=np.shape(data)[1] # number of data samples per trial
    k=np.shape(data)[2] # number of trials
    time=N/sampling_rate # timelength of data in seconds
    
#    The update coefficient Uc (0<Uc<1)
#    controls the adaptation speed of time-varying MVAR parameters             
    Uc=0.02
    
    Impone=np.eye(m*p)
    Imone=np.eye(m)
    Ikone=np.eye(k)
    
    Awave=np.zeros((m*p,m,N))
    Hp=np.zeros((k,m*p,N))
    
    P=np.tile(np.eye(m*p), [1,1,N])
    W=np.tile(np.eye(m), [1,1,N])
    E=np.zeros((k,m,N))
    X=np.zeros((k,k,N))
    KG=np.zeros((m*p,k,N))
    V=np.zeros((m*p,m*p,N))
    
    for n in range(p,N):
        for i in range(0,p):
            O=np.transpose(np.squeeze(data[:][n-i][:]))
            Hp[:][(i-1)*m+1:(i*m)+1][n]=O
            
        dataE=np.squeeze(data[:][n][:])
        E[:][:][n]=np.transpose(dataE)-Hp[:][:][n]*Awave[:][:][n-1]
        W[:][:][n]=(1-Uc)*W[:][:][n-1]+Uc*(np.transpose(E[:][:][n])*E[:][:][n])/(k-1)
        X[:][:][n]=np.invert(Hp[:][:][n]*P[:][:][n-1]*np.transpose(Hp[:][:][n])+ np.trace(W[:][:][n])*Ikone)
        KG[:][:][n]=P[:][:][n-1]*np.transpose(Hp[:][:][n])*X[:][:][n]
        Awave[:][:][n]=Awave[:][:][n-1]+KG[:][:][n]*E[:][:][n]
        V[:][:][n]=((Uc*np.trace((Impone-KG[:][:][n]*Hp[:][:][n])*P[:][:][n-1]))/(m*m*p))*Impone
        P[:][:][n]=(Impone-KG[:][:][n]*Hp[:][:][n])*P[:][:][n-1]+V[:][:][n]
        
    ApAR=np.zeros((m,m,p,p))
    for i in range(p,N):
        for q in range(0,p):
            ApAR[:][:][q][i]=np.transpose(Awave[(m*(q-1)+1):m*q][:][i])
            
            
    return {'ApAR':ApAR,  'E':E, 'nchannels':m, 'Order':p, 'nsamples':N, 'ntrials':k, 'sampling_rate':sampling_rate, 'time':time }
        