# scripts to compare fit files from fitacf and fitlomb

# generate rawacf file with known targets
# process file with fitacf and fitlomb
# compare output..

import os
import sys
import pdb
import datetime
import h5py

from multiprocessing import Manager
from rawacf_generator import rawacf_record
from superdarn_tools import *

sys.path.append('../SuperDARN_FitLOMB')
from pydarncuda_fitlomb import generate_fitlomb 
SANDBOX = '/home/kleinjt/repos/SuperDARN_pydmap_write/sandbox'
ACF_NAME = '20150101.0000.00.tst'
RAWACF_EXT = '.rawacf'
FITACF_EXT = '.fitacf'
RADAR = 'tst'
RECORDTIME = datetime.datetime(2015, 1, 1, 0, 0)
DEF_NOISE = 1

class target(object):
    # defines a target a range rangekm with velocity and width in meters per second 
    def __init__(self, rangegate = 0, velocity = 0, width = 0, power = 0, v_e = 0, w_e = 0, nlag = 0):
        self.rangegate = rangegate
        self.velocity = velocity
        self.width = width
        self.power = power
        self.v_e = v_e
        self.w_e = w_e
        self.nlag = nlag


# synthesizes a rawacf file with velocity and spectral width at gates rgates 
def generate_rawacf(rawacfname, targets = [], noise = DEF_NOISE):
    rawacf = rawacf_record(filename = rawacfname)

    for target in targets:
        rawacf.addTarget(target.rangegate, target.power, target.velocity, target.width)
    
    rawacf.generateScatter()
    rawacf.applyNoise(noise)
    rawacf.setTime(RECORDTIME) 
    rawacf.write()
    rawacf.close()

# generates a fitacf file fitacfname from rawacf file rawacfname
def generate_fitacf(rawacfname, fitacfname):
    cmd = 'make_fit -new {} > {}'.format(rawacfname, fitacfname)
    os.system(cmd)

# generates a fitlomb file fitacfname from rawacf file rawacfname
def generate_fitlomb(rawacfname):
    os.environ['DAVIT_LOCALDIR'] = SANDBOX
    os.environ['DAVIT_DIRFORMAT'] = '%(dirtree)s/'
    radar = RADAR
    manager = Manager()
    lock = manager.Lock()

    stime = RECORDTIME
    etime = None
    
    record = ((stime, etime, radar, lock))
    generate_fitlomb(record)

# parse fitacf, returns a list of targets at ranges rgates
def parse_fitacf(fitacfname, rgates):
    dumpstr = get_dmapdumpstring()
    scandata = parse_dmapdumpstring(dumpstr)
    targets = []
    for rgate in rgates:
        if rgate in scandata['slist']:
            sidx = np.nonzero(scandata['slist'] == rgate)[0][0]
            v = scandata['v'][sidx]
            w_l = scandata['w_l'][sidx]
            p_l = scandata['p_l'][sidx]
            v_e = scandata['v_e'][sidx]
            w_l_e = scandata['w_l_e'][sidx]
            nlag = scandata['nlag'][sidx]

            targets.append(target(rangegate = rgate, velocity = v, width = w_l, power = p_l, v_e = v_e, w_e = w_l_e, nlag = nlag))

        else:
            print('scatter not found on rage gate {}'.format(rgate))
            targets.append(target(rangegate = rgate, nlag = 0))

    return targets

# parse fitlomb, return a list of targets at ranges rgates
def parse_fitlomb(fitlombname, rgates):
    h5f = h5py.File(fitlombname, 'r')
    beamgrp = h5f[h5f.keys()[1]]
    targets = []

    v = h5f[beamgrp]['v']
    v_e = h5f[beamgrp]['v_e']
    p_l = h5f[beamgrp]['p_l']
    w_l = h5f[beamgrp]['w_l']
    w_l_e = h5f[beamgrp]['w_l_e']
    nlag = h5f[beamgrp]['nlag']

    for rgate in rgates:
        targets.append(target(rangegate = rgate, velocity = v[rgate], width = w_l[rgate], power = p_l[rgate], v_e = v_e[rgate], w_e = w_l_e[rgate], nlag = nlag))

    return targets

def main():
    rawacf_name = SANDBOX + '/' + ACF_NAME + RAWACF_EXT
    fitacf_name = SANDBOX + '/' + ACF_NAME + FITACF_EXT
    synthetic_targets = []
    t = target(rangegate = 5, velocity = 500, width = 200, power = 1)
    synthetic_targets.append(t)


    gates = [t.rangegate for t in synthetic_targets]
    
    generate_rawacf(rawacf_name, targets = synthetic_targets, noise = .1)
    generate_fitacf(rawacf_name, fitacf_name)
    generate_fitlomb(rawacf_name)
    fitacf_targets = parse_fitacf(fitacf_name, gates)
    #fitlomb_targets = parse_fitlomb(fitlomb_name, gates)
        
    pdb.set_trace()
    



if __name__ == '__main__':
    main()
