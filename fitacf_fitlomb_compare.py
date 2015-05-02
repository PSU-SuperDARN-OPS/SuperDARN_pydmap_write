# scripts to compare fit files from fitacf and fitlomb

# generate rawacf file with known targets
# process file with fitacf and fitlomb
# compare output..

import os
import sys
import pdb
import datetime
import h5py
import pydarn.sdio as sdio

from multiprocessing import Manager
from rawacf_generator import rawacf_record
from superdarn_tools import *

sys.path.append('../SuperDARN_FitLOMB')
import lagstate
from cuda_bayes import BayesGPU
from pydarncuda_fitlomb import CULombFit

SANDBOX = '/home/radar/repos/SuperDARN_pydmap_write/sandbox'
ACF_NAME = '20150101.0000.00.tst'
RAWACF_EXT = '.rawacf'
FITACF_EXT = '.fitacf'
HDF5_EXT = '.fitlomb.hdf5'

RADAR = 'tst'
RECORDTIME = datetime.datetime(2015, 1, 1, 0, 0)
DEF_NOISE = .1

C = 3e8
MAX_TFREQ = 16e6
LOMB_PASSES = 1
NFREQS = 256
NALFS = 256
FWHM_TO_SIGMA = 2.355 # conversion of fwhm to std deviation, assuming gaussian
MAX_V = 2000 # m/s, max velocity (doppler shift) to include in lomb
MAX_W = 1200 # m/s, max spectral width to include in lomb 


# synthesizes a rawacf file with velocity and spectral width at gates rgates 
def generate_rawacf(rawacfname, targets = [], noise = DEF_NOISE):
    rawacf = rawacf_record(filename = rawacfname)

    for target in targets:
        rawacf.addTarget(target)
    
    rawacf.generateScatter()
    rawacf.applyNoise(noise)
    rawacf.calcPwr0()
    rawacf.setTime(RECORDTIME) 
    rawacf.write()
    rawacf.close()

# generates a fitacf file fitacfname from rawacf file rawacfname
def generate_fitacf(rawacfname, fitacfname):
    cmd = 'make_fit -new {} > {}'.format(rawacfname, fitacfname)
    os.system(cmd)

# generates a fitlomb file fitacfname from rawacf file rawacfname
def generate_fitlomb(rawacfname, fitlombname):
    os.environ['DAVIT_LOCALDIR'] = SANDBOX
    os.environ['DAVIT_DIRFORMAT'] = '%(dirtree)s/'
    radar = RADAR
    manager = Manager()
    lock = manager.Lock()
    
    stime = RECORDTIME
    etime = None
    
    hdf5file = h5py.File(fitlombname, 'w')

    myPtr = sdio.radDataOpen(stime,radar,eTime=etime,channel=None,bmnum=None,cp=None,fileType='rawacf',filtered=False, src='local', noCache = True)
    drec = sdio.radDataReadRec(myPtr)

    amax = np.ceil((np.pi * 2 * MAX_TFREQ * MAX_W) / C)
    fmax = np.ceil(MAX_V * 2 * MAX_TFREQ / C)
    freqs = np.linspace(-fmax,fmax, NFREQS)
    alfs = np.linspace(0, amax, NALFS)

    fit = CULombFit(drec)
    gpu_lambda = BayesGPU(fit.lags, freqs, alfs, fit.nrang, LAMBDA_FIT)
    txlag_cache = lagstate.good_lags_txsamples(fit)
    fit.SetBadlags(txlag_cache = txlag_cache)
    fit.CudaProcessPulse(gpu_lambda)
    fit.CudaCopyPeaks(gpu_lambda)
    fit.WriteLSSFit(hdf5file)
    hdf5file.close()

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
            print('scatter not found on range gate {}'.format(rgate))
            targets.append(target(rangegate = rgate, nlag = 0))

    return targets

# parse fitlomb, return a list of targets at ranges rgates
def parse_fitlomb(fitlombname, rgates):
    h5f = h5py.File(fitlombname, 'r')
    beamgrp = h5f[h5f.keys()[0]]
    targets = []

    v = beamgrp['v'][...]
    v_e = beamgrp['v_e'][...]
    p_l =  beamgrp['p_l'][...]
    w_l =  beamgrp['w_l'][...]
    w_l_e = beamgrp['w_l_e'][...]
    nlag = beamgrp['nlag'][...]
    for rgate in rgates:
        targets.append(target(rangegate = rgate, velocity = v[rgate], width = w_l[rgate], power = p_l[rgate], v_e = v_e[rgate], w_e = w_l_e[rgate], nlag = nlag[rgate]))

    return targets

def test_fitacf():
    rawacf_name = SANDBOX + '/' + ACF_NAME + RAWACF_EXT
    fitacf_name = SANDBOX + '/' + ACF_NAME + FITACF_EXT
    synthetic_targets = []
    t = target(rangegate = 5, velocity = 500, width = 200, power = 1)
    synthetic_targets.append(t)

    gates = [t.rangegate for t in synthetic_targets]
    generate_rawacf(rawacf_name, targets = synthetic_targets, noise = .1)
    generate_fitacf(rawacf_name, fitacf_name)
    fitacf_targets = parse_fitacf(fitacf_name, gates)
    
    print 'synthetic targets:'
    for t in synthetic_targets:
        print t
    print 'fitacf fitted targets:'
    for t in fitacf_targets:
        print t


def main():
    print 'finished importing...'
    rawacf_name = SANDBOX + '/' + ACF_NAME + RAWACF_EXT
    fitacf_name = SANDBOX + '/' + ACF_NAME + FITACF_EXT
    fitlomb_name = SANDBOX + '/' + ACF_NAME + HDF5_EXT

    print 'creating targets...'
    synthetic_targets = []
    synthetic_targets.append(target(rangegate = 5, velocity = 500, width = 200, power = 1))
    synthetic_targets.append(target(rangegate = 6, velocity = 300, width = 100, power = 1))
    synthetic_targets.append(target(rangegate = 7, velocity = 000, width = 000, power = 1))


    gates = [t.rangegate for t in synthetic_targets]
    
    print 'creating rawacf...'
    generate_rawacf(rawacf_name, targets = synthetic_targets, noise = .1)

    print 'creating fitacf...'
    generate_fitacf(rawacf_name, fitacf_name)
    
    print 'creating fitlomb...'
    generate_fitlomb(rawacf_name, fitlomb_name)

    fitacf_targets = parse_fitacf(fitacf_name, gates)
    fitlomb_targets = parse_fitlomb(fitlomb_name, gates)

    print 'synthetic targets:'
    for t in synthetic_targets:
        print t

    print 'fitacf fitted targets:'
    for t in fitacf_targets:
        print t

    print 'fitlomb fitted targets:'
    for t in fitlomb_targets:
        print t

if __name__ == '__main__':
    main()
