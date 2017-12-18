# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 13:14:23 2017

@author: Administrator
"""
import numpy as np
from matplotlib import pyplot

def fcon_value(T=None,n=None,rf=None,rq=None,s0=None,strike=None,sigma=None,C=None,SB=None):
    '''
    #T expiration date
    #n  total number of periods
    #stike exercise  price
    #sigma voltaility of 
    #s0 present price 
    #rf risk-free rate
    #rq credit rate
    #B book value
    #C redemption price
    #SB sell back price
    '''
    deltat=T/n
    u=np.exp(sigma*deltat)
    d=1/u
    p=(np.exp(rf*np.sqrt(deltat))-d)/(u-d)
    q=1-p
    B=100
    k0=B/strike
    bound=strike*0.9
    
    s=np.zeros([n,n+1])
    v=np.zeros([n,n+1])
    D=np.zeros([n,n+1])
    r=np.zeros([n,n+1])
    #与原matlab代码不同修改了0轴维度由n+1变为n
    k=np.zeros([n,n+1])
    
    #构造s
    for j in np.arange(0,n).reshape(-1):
        for i in np.arange(0,j+2).reshape(-1):
            s[j,i]=s0*(u**(j+1-i))*(d**(i))
    
    k[:,0]=k0#此处与原matlab代码有改动
    k[0,1]=k0
    
    #构造k    
    for j in np.arange(1,n).reshape(-1):
        for i in np.arange(1,j+2).reshape(-1):
            if k[j-1,i-1] == k0:
                if s[j,i]<bound:
                    k[j,i]=B/s[j,i]
                else:
                    k[j,i]=k0
            else:
                k[j,i]=k[j-1,i-1]
    
    #构造初始的第n层,第n行表示第n期的结果
    for i in np.arange(0,n+1):
        v[n-1,i]=np.max([k[n-1,i]*s[n-1,i],B])
        if v[n-1,i] > B:
            r[n-1,i]=rf
        else:
            r[n-1,i]=rq
    
    #够造r
    for j in np.arange(n-2,-1,-1).reshape(-1):
        for i in np.arange(0,j+2).reshape(-1):
            r[j,i]=p*r[j+1,i]+q*r[j+1,i+1]
            if k[j,i]*s[j,i]>B:
                r[j,i]=rf
    
    #第n-1层
    for i in np.arange(0,n):
        v[n-2,i]=(p*v[n-1,i]+q*v[n-1,i+1])*np.exp(-r[n-2,i]*deltat)
        D[n-2,i]=max(k[n-2,i]*s[n-2,i],min(v[n-2,i],C),SB)
    
    for j in np.arange(n-3,-1,-1).reshape(-1):
        for i in np.arange(0,j+2):
            v[j,i]=(p*D[j+1,i]+q*D[j+1,i+1])*np.exp(-r[j,i]*deltat)
            D[j,i]=max(k[j,i]*s[j,i],min(v[j,i],C),SB)
    
    r0=p*r[0,0]+q*r[0,1]
    value=(p*D[0,0]+q*D[0,1])*np.exp(-r0*deltat)
    
    pyplot.hist(D[n-2,:-1],20,normed=1)
    pyplot.show()
    
    print('the price is %f'%(value))
    return value
    
if __name__ == '__main__':
    p=fcon_value(T=5,n=10,rf=0.03,rq=0.05,s0=50,strike=50,sigma=0.3,C=110,SB=85)
        
    
    
    