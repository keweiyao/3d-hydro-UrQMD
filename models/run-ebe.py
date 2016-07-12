import subprocess as sp
import numpy as np
import sys
import h5py as h5
import os

mult_low = 150.0
mult_high = 300.0

# arguments: event number, will be modified to the random seed of the event
event_number = int(sys.argv[1])

# Initial condition: TRENTo model
os.chdir("./trento")

np.random.seed(event_number)
random_seed = np.random.randint(0, 1000000)
print "try random seed = ", random_seed
cmd = "./trento -c trento-config.txt --random-seed %d"%random_seed
sp.call(cmd, shell=True)
f = h5.File("./ic.hdf5")
attr = f[f.keys()[0]].attrs
mult = attr['mult']

while mult < mult_low or mult > mult_high:
    cmd = "rm ./ic.hdf5"
    sp.call(cmd, shell=True)

    random_seed = np.random.randint(0, 1000000)
    print "try random seed = ", random_seed
    cmd = "./trento -c trento-config.txt --random-seed %d"%random_seed
    sp.call(cmd, shell=True)
    f = h5.File("./ic.hdf5")
    attr = f[f.keys()[0]].attrs
    mult = attr['mult']
print "Event selection succeed, mult = %f, b = %f", mult, attr['b']
os.chdir("..")




# Determine how many oversamples are needed by Number of participants.
N_samples_basic = 5 #number of basic samples
N_samples_max = 400 #number of maximum samples
N_part_max = 208*2 #maximum number of participants of Pb+Pb
f = h5.File("./trento/ic.hdf5")
attr = f[f.keys()[0]].attrs
npart = attr['npart']
print "Npart = ", npart
Ns = np.min([int(N_samples_basic*N_part_max/npart), N_samples_max])
print "Oversamples = ", Ns
cmd = 'sed -i "s/number_of_events.*/number_of_events %d/g" ./escape/escape-config.txt'%Ns
sp.call(cmd, shell=True)

# Move ic file to hydro and run hydrodynamics
cmd = "mv -v ./trento/ic.hdf5 ./vhlle"
sp.call(cmd, shell=True)
os.chdir("./vhlle")
cmd = "./vhlle vhlle-config.txt"
sp.call(cmd, shell=True)
os.chdir("..")
cmd = "mv -v ./vhlle/vhlle_result_dir/freezeout.dat ./escape"
sp.call(cmd, shell=True)

# Run afterburner with adaptive number of oversampling
os.chdir("./escape")
cmd = "./escape events %d escape-config.txt > afterburner%d.dat"%(event_number, event_number)
sp.call(cmd, shell=True)
os.chdir("..")

# Collect data
os.chdir("./collector")
cmd = "cp ../escape/afterburner%d.dat ./plist.dat"%event_number
sp.call(cmd, shell=True)
cmd = "./collector stat-config.txt"
sp.call(cmd, shell=True)

# calculate some observables
cmd = "mkdir data-folder && cp data.hdf5 data-folder"
sp.call(cmd, shell=True)
cmd = "python calc_obs.py %d"%event_number
sp.call(cmd, shell=True)
os.chdir("..")

# prepare results
cmd = "mv -v collector/data.hdf5 ./plist%d.hdf5 && mv -v collector/Obs.dat ./Obs%d.dat"%(event_number, event_number)
sp.call(cmd, shell=True)

# done
print "Event ", event_number, "done!"
