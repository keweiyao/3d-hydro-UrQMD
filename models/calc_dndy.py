import numpy as np
import h5py
import os
import sys
import scipy.interpolate as interp

def extract_dNdy(filename):
	f = h5py.File(filename, 'r')
	dNdy = []
	ET = []
	if len(f.keys()) < 1:
		return False, np.zeros(25)
	for sample in f.keys():
		eta = f[sample]['eta'].value		
		ch = f[sample]['ch'].value

		eta = eta[ch!= 0]
		dndy, xbins = np.histogram(eta, bins = 25, range = [-6, 6])
		dNdy.append(dndy)
	dNdy = np.array(dNdy)
	print ((xbins[1:] + xbins[:-1])*0.5)
	return True, dNdy

def extract_dsdy(filename):
	f = h5py.File(filename, 'r')
	data = []
	N_events = len(f.keys())
	cen = (np.array([0., 5., 10., 20., 30., 40., 50., 60., 70., 80])/100.*N_events).astype(int)
	for k in f.keys():
		dsdy = f[k].value[0,0,:]
		ymax = 10.
		Ny = dsdy.shape[0]
		y = np.linspace(-ymax, ymax, Ny)
		dy = y[1] - y[0]
		#dxdy = 0.1**2
		#dsdy = np.array([np.sum(ds[:,:,i]) for i in range(Ny)])*dxdy
	
		yout = np.linspace(-6, 6, 26)
		yout = 0.5*(yout[1:]+yout[:-1])
		F = interp.interp1d(y-dy*0.5, dsdy)
		output = F(yout)
		data.append(output)
	sorteddata = np.array(sorted(data, key=lambda x: np.sum(x[11:14]), reverse=True))
	obs = [np.mean(sorteddata[cen[i]:cen[i+1],:], axis=0) for i in range(9)]
	return np.array(obs)
	
fout = h5py.File("dndy.hdf5", 'w')
#fin_ic = sys.argv[1]
fin_hybrid = sys.argv[1]

status, dndy = extract_dNdy(fin_hybrid)
#dsdy = extract_dsdy(fin_ic)

#fout.create_dataset('dsdy', data=dsdy)
fout.create_dataset('dNdy', data=dndy)
fout.create_dataset('status', data=status)
fout.close()
