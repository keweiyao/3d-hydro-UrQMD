import numpy as np
import h5py
import sys

filename = sys.argv[1]
count = 0
fout = h5py.File('final-particles.hdf5', 'w')

class sample():
    __slots__ = ["sample_id", "pt", "eta", "E", "phi", "pid", "ch"]
    def __init__(self, sid):
        self.sample_id = sid
        self.pt = []
        self.eta = []
        self.E = []
        self.phi = []
        self.pid = []
        self.ch = []

sample_list = []
with open(filename) as f:
    for l in f:
        nl = l.split()
        if nl[0] == "#event":
            this_sample = sample(int(nl[1]))
            sample_list.append(this_sample)
            count = 0
        if nl[0] == "#cols:":
            count += 1
        if count == 2 and len(nl) == 13:
            px = float(nl[4])
            py = float(nl[5])
            pz = float(nl[6])
            E = float(nl[7])
            pid = int(nl[8])
            ch = int(nl[10])
            pt = (px**2 + py**2)**0.5
            phi = np.arctan2(py, px)
            pabs = (pt**2 + pz**2)**0.5
            eta = 0.5*np.log((pabs+pz)/(pabs-pz))
            
            sample_list[-1].pt.append(pt)
            sample_list[-1].E.append(E)
            sample_list[-1].eta.append(eta)
            sample_list[-1].phi.append(phi)
            sample_list[-1].pid.append(pid)
            sample_list[-1].ch.append(ch)


for sample in sample_list:
    group = fout.create_group("sample-%d"%sample.sample_id)
    print (group.name)
    group.create_dataset("pt", data=sample.pt)
    group.create_dataset("eta", data=sample.eta)
    group.create_dataset("E", data=sample.E)
    group.create_dataset("phi", data=sample.phi)
    group.create_dataset("pid", data=sample.pid)
    group.create_dataset("ch", data=sample.ch)
fout.close()
