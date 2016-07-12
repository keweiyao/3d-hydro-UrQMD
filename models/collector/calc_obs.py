import numpy as np
import h5py as h5
import os
import sys

def calc_v1(ds, cut):
    pt = ds['pt']
    psin = ds['psin']
    pcos = ds['pcos']
    eta = ds['eta']
    s = psin[(cut[0] < pt) & (pt < cut[1]) & (eta < cut[3]) & (eta>cut[2])]
    c = pcos[(cut[0] < pt) & (pt < cut[1]) & (eta < cut[3]) & (eta>cut[2])]  
    w = np.size(s)
    if w > 0:
        return w, np.sum(c)/w, np.sum(s)/w
    else:
        return w, 0.0, 0.0
def calc_cn2(ds, cut):
    pt = ds['pt']
    psin = ds['psin']
    pcos = ds['pcos']
    eta = ds['eta']
    s = psin[(cut[0] < pt) & (pt < cut[1]) & (eta < cut[3]) & (eta>cut[2])]
    c = pcos[(cut[0] < pt) & (pt < cut[1]) & (eta < cut[3]) & (eta>cut[2])]                                                                           
    c2 = c*c;          s2 = s*s;           cs = c*s;                                                               
    c3 = c2*c;         c2s = c2*s;         cs2 = c*s2;        s3 = s2*s;                                           
    c4 = c3*c;         c3s = c3*s;         c2s2 = c2*s2;      cs3 = c*s3;       s4 = s3*s;                         
    c5 = c4*c;         c4s = c4*s;         c3s2 = c3*s2;      c2s3 = c2*s3;     cs4 = c*s4;       s5 = s4*s;
    c6 = c5*c;         c5s = c5*s;         c4s2 = c4*s2;      c3s3 = c3*s3;    c2s4 = c2*s4;     cs5 = c*s5;   s6 = s5*s;
    Q1 = np.sum(c+1j*s)                  
    Q2 = np.sum(c2 - s2 + 1j*2.0*cs)                
    Q3 = np.sum(c3 - 3.0*cs2 + 1j*( 3.0*c2s - s3 ))                        
    Q4 = np.sum(c4 + s4 - 6.0*c2s2 + 1j*4.0*(c3s - cs3))                    
    Q5 = np.sum(c5 - 10.0*c3s2 + 5.0*cs4 + 1j*(s5+5.0*c4s - 10.0*c2s3) )
    Q6 = np.sum(c6 - s6 + 15.0*c2s4 - 15.0*c4s2 + 1j*(6.0*c5s - 20.0*c3s3 + 6.0*cs5))       
    M = np.size(s)
    
    w2 = M*(M-1)
    c22 = 0.0
    c32 = 0.0
    c42 = 0.0
    c52 = 0.0
    c62 = 0.0
    if w2 != 0:
        c22 = (np.abs(Q2)**2 - M)/w2
        c32 = (np.abs(Q3)**2 - M)/w2
        c42 = (np.abs(Q4)**2 - M)/w2
        c52 = (np.abs(Q5)**2 - M)/w2
        c62 = (np.abs(Q6)**2 - M)/w2
    return w2, c22, c32, c42, c52, c62

def calc_dndy(ds, cut):
    pt = ds['pt']
    psin = ds['psin']
    eta = ds['eta']
    p = psin[(cut[0] < pt) & (pt < cut[1]) & (eta < cut[3]) & (eta>cut[2])]
    return [np.size(p)*1.0/(cut[3]-cut[2])]

def calc_cen(ds, cut):
    pt = ds['pt']
    psin = ds['psin']
    eta = ds['eta']
    p = psin[(cut[0] < pt) & (pt < cut[1]) & (eta < cut[3]) & (eta>cut[2])]
    return np.size(p)

event_number = int(sys.argv[1])
deta_cn2 = 1.6
deta_dndy = 1.0
deta_v1 = 0.5
fo = open('Obs.dat','w')
dict = {}
if event_number >= 0:
    f = h5.File("data-folder/data.hdf5", 'r')
    k = f.keys()
    # calculate and append two-particle correlation
    cn2 = []
    v1 = []
    dndy = []
    cen_key = []
    for c in range(9):
        cut = [0.2, 5.0, (2.*c-9)/2.*deta_cn2, (2.*c-7)/2.*deta_cn2]
        result = []
        for ik in k:
            ds = f[ik].value
            result.append(calc_cn2(ds, cut))
        result = np.array(result)
        w = result[:,0]
        sample_avg = np.average(result[:, 1:], axis = 0, weights = w) 
        cn2.append(np.concatenate((cut, [np.mean(w)], sample_avg)))
#    print "cn2, ", cn2
    cn2 = np.array(cn2)
    dict.update({'cn2' : cn2})
    
    for c in range(9):
        cut = [0.15, 5.0, (2.*c-9)/2.*deta_v1, (2.*c-7)/2.*deta_v1]
        result = []
        for ik in k:
            ds = f[ik].value
            result.append(calc_v1(ds, cut))
        result = np.array(result)
        w = result[:,0]
        sample_avg = np.average(result[:, 1:], axis = 0, weights = w)
        v1.append(np.concatenate((cut, [np.mean(w)], sample_avg)))
    v1 = np.array(v1)
#    print "v1, ", v1
    dict.update({'v1' : v1})

    for c in range(13):
        cut = [0.0, 100.0, (2.*c-13)/2.*deta_dndy, (2.*c-11)/2.*deta_dndy]
        result = []
        for ik in k:
            ds = f[ik].value
            result.append(calc_dndy(ds, cut))
        result = np.array(result)
        sample_avg = np.mean(result, axis = 0)
        dndy.append(np.concatenate((cut, sample_avg)))
    dndy = np.array(dndy)
#    print "dndy, ", dndy
    dict.update({'dndy' : dndy})
    
    cut = [0, 100.0, -0.8, 0.8]
    result = []
    for ik in k:
        ds = f[ik].value
        result.append(calc_cen(ds, cut))
    result = np.array(result)
    sample_avg = np.mean(result)
    cen_key = np.concatenate((cut, [sample_avg]))
#    print "cen_key, ", cen_key
    dict.update({'cen' : cen_key})

np.savez(fo, **dict)
fo.close()
