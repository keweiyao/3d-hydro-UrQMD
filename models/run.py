import subprocess as sp
import numpy as np
import sys
import h5py as h5
import os

def cmd(acmd):
	sp.call(acmd, shell=True)
def cd(dir):
	os.chdir(dir)

oversample_scheme = 1
Basic_oversamples = 50
Max_oversamples = 500
config_file = sys.argv[1]
ic_file = sys.argv[2]
event_number = int(sys.argv[3])

np.random.seed(event_number)
random_seed = np.random.randint(1)
# arguments: event number
#------------Separate the large configuration file into small ones---------
cmd(" awk -v RS= '{print > NR}' ./%s "%config_file)
cmd(" sed -e '1,1d' 1 > ./trento3d/trento-config.txt ")
cmd(" sed -e '1,1d' 2 >  ./vhlle/vhlle-config.txt ")
cmd(" sed -e '1,1d' 3 >  ./escape/escape-config.txt ")

cmd("sed -i -e 's/=/ /g' ./vhlle/vhlle-config.txt") 
cmd("sed -i -e 's/=/ /g' ./escape/escape-config.txt")

#------------Initial condition: TRENTo model-----------------------
"""
cd("./trento3d")
cmd("./trento3d Pb Pb 1000 -c trento-config.txt --random-seed %d"%(random_seed))
f = h5.File("./ic.hdf5")
# averageing entropy density
mult = []
Npart = []
Ncen = (np.array([0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.])/100.*1000.).astype(int)
ds = []
key = f.keys()
for k in f.keys():
    attr = f[k].attrs
    mult.append(attr['mult']); Npart.append(attr['npart'])
mult = np.array(mult)
Npart = np.array(Npart)
ds = np.array(ds)
index = np.argsort(mult)[::-1]
#f.close()

ofilename = "%d-avg-ic.hdf5"%event_number
f2 = h5.File(ofilename, 'w')
for i in range(10):
    print np.mean(mult[index[Ncen[i]:Ncen[i+1]]]), np.mean(Npart[index[Ncen[i]:Ncen[i+1]]]) 
    mM, mN = np.mean(mult[index[Ncen[i]:Ncen[i+1]]]), np.mean(Npart[index[Ncen[i]:Ncen[i+1]]])
    subevents = np.mean(np.array([ f[key[j]].value for j in index[Ncen[i]:Ncen[i+1]] ]), axis=0)
    dataset = f2.create_dataset('cen-avg-%d'%i, data=subevents)
    dataset.attrs.create("Npart", data=mN)
    dataset.attrs.create("Mult", data=mM)
f.close()
f2.close()
#print ("IC: \n  mean mult = %f \n mean npart = %f"%(mM, mN) )
cd("..")

cmd("mv -v ./trento3d/%s ./"%ofilename)

"""

#cmd("mv -v ./%s ./vhlle/ic.hdf5"%ic_file)
f = h5.File(ic_file)
ds = f['cen-avg-%d'%event_number].value*1.2
Npart = f['cen-avg-%d'%event_number].attrs['Npart']

f.create_dataset('event_0', data=ds)
f.close()
#------------Determine over samples--------------------------------
if oversample_scheme == 0: # fixed number of oversamples
	Ns = Basic_oversamples
if oversample_scheme == 1: # fixed total charged particles
	Ns = int(Basic_oversamples*408.0/Npart)
	Ns = np.min([Ns, Max_oversamples])
print ("Oversamples = %d"%Ns)
cmd('sed -i "s/number_of_events.*/number_of_events %d/g" ./escape/escape-config.txt'%Ns)
cmd("mv -v ./%s ./vhlle/ic.hdf5"%ic_file)

#----------------------vhlle---------------------------------------
# Move ic file to hydro and run hydrodynamics

#cmd("mv -v ./trento3d/ic.hdf5 ./vhlle")
cd("./vhlle")
cmd("./vhlle vhlle-config.txt")
cd("..")

#----------------------Escape5-------------------------------------
cmd("mv -v ./vhlle/vhlle_result_dir/freezeout.dat ./escape")
cd("./escape")
cmd("./escape events %d escape-config.txt > afterburner.dat"%(event_number))
os.chdir("..")

#---------------------Data-preposessing----------------------------
cmd("python preprocessing.py  escape/afterburner.dat ")

#---------------------Optional-Precalculate-dndy-------------------
cmd("python calc_dndy.py final-particles.hdf5")

#-----------------prep-data-for-transfer---------------------------
#cmd("mv final-particles.hdf5  final-particles-%d.hdf5 "%event_number)
cmd("mv dndy.hdf5  dndy-%d.hdf5 "%(event_number) )

#  finished
print ("Event done!")
