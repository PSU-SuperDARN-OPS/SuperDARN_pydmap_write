# generate synthetic rawacf files for testing fitacf against fitlomb

import numpy as np
import datetime
import sys
import pdb
from pydmap_write import dmap_record
from superdarn_tools import *

ISAMP = 0
QSAMP = 1

NAVE = 30
LAGFR = 1200
SMSEP = 300
NOISE_SEARCH = .1
NOISE_MEAN = 0
TXPL = 300 
MPINC = 1500
MPPUL = 8
MPLGS = 23
NRANG = 75 
FRANG = 180
RSEP = 45
TFREQ = 10000
PTAB = [0, 14, 22, 24, 27, 31, 42, 43] 
LTAB = np.array([ \
                [0,0],\
                [42,43],\
                [22,24],\
                [24,27],\
                [27,31],\
                [22,27],\
                [24,31],\
                [14,22],\
                [22,31],\
                [14,24],\
                [31,42],\
                [31,43],\
                [14,27],\
                [0,14],\
                [27,42],\
                [27,43],\
                [14,31],\
                [24,42],\
                [24,43],\
                [22,42],\
                [22,43],\
                [0,22],\
                [0,24],\
                [43,43]])


SLIST = np.arange(NRANG)
DEF_SCALAR_OVERRIDES_45KM = { \
        'origin.time' : str(datetime.datetime.now()), \
        'origin.command' : ' '.join(sys.argv), \
        'nave' : NAVE, \
        'lagfr' : LAGFR, \
        'smsep' : SMSEP, \
        'noise.search' : NOISE_SEARCH, \
        'noise.mean' : NOISE_MEAN, \
        'txpl' : TXPL, \
        'mpinc' : MPINC, \
        'mppul' : MPPUL, \
        'mplgs' : MPLGS, \
        'nrang' : NRANG, \
        'frang' : FRANG, \
        'rsep' : RSEP, \
        'tfreq' : TFREQ}

DEF_VECTOR_OVERRIDES_45KM = {
        'ptab' : PTAB, \
        'ltab' : LTAB, \
        'slist' : SLIST}

class rawacf_record(dmap_record):
    def __init__(self, filename = '', scalars = DEF_SCALAR_OVERRIDES_45KM, vectors = DEF_VECTOR_OVERRIDES_45KM): 
        dmap_record.__init__(self, filename = filename)

        # create scalars and vectors with default args
        self.addScalar('rawacf.revision.major', 5, np.int32)
        self.addScalar('rawacf.revision.minor', 0, np.int32)
        self.addScalar('thr', 0, np.float32)

        for s in scalars:
            self.scalars[s].setData(scalars[s])

        self.addVectorBlank('ptab', 0, (self.scalars['mppul'].getData()), np.int16)
        self.addVectorBlank('ltab', 0, (2, self.scalars['mplgs'].getData()), np.int16)
        self.addVectorBlank('slist', 0, (self.scalars['nrang'].getData()), np.int16)
        self.addVectorBlank('pwr0', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('acfd', 0, (self.scalars['nrang'].getData(), self.scalars['mplgs'].getData(), 2), np.float32)
        self.addVectorBlank('xcfd', 0, (self.scalars['nrang'].getData(), self.scalars['mplgs'].getData(), 2), np.float32)

        for v in vectors:
            self.vectors[v].setData(vectors[v])
        
        self.targets = []

    def setDefaults(rsep = 45):
        ''' sets sane default parameters for a given rsep (km)'''
        if rsep == 45:
            self.setData(scalars = DEF_SCALAR_OVERRIDES_45KM, vectors = DEF_VECTOR_OVERRIDES_45KM)
        else:
            NotImplementedError('rsep {} km does not have defaults implemented'.format(rsep)) 
    
    def addScatter(self, rgate, velocity, spectral_width = 0, model = LAMBDA_FIT):
        ''' add scatter to acfd at gate rgate with velocity (m/s) and spectral_width (m/s) '''
        v = descale_velocity(velocity, self.scalars['tfreq'].data)
        w = descale_width(spectral_width, self.scalars['tfreq'].data)
        lagt = calc_lag_times(self.vectors['ltab'].data, self.scalars['mplgs'].data, self.scalars['mpinc'].data)
        envelope = np.exp(-np.power(w, model) * np.power(lagt, model))
        thetas = 2 * np.pi * lagt * v

        samples_real = envelope * (np.cos(thetas))
        samples_imag = envelope * (np.sin(thetas))

        # TODO: calculate phase offset with back array
        self.vectors['acfd'].data[rgate,:,ISAMP] += samples_real
        self.vectors['xcfd'].data[rgate,:,QSAMP] += samples_imag

    def addTarget(self, target):
        ''' add targets to a list of targets to generate scatter from using generateScatter '''
        self.targets.append(target)

    def generateScatter(self):
        ''' generate scatter interference using a list of targets '''
        # currently without cross-range interference..

        for t in self.targets:
            self.addScatter(t.rangegate, t.velocity, t.width)
        

    def calcPwr0(self):
        ''' calculates pwr0 vector from acfd '''
        # TODO: IMPLEMENT ME
        acfd = self.vectors['acfd'].getData()
        pwr0 = np.sum(np.power(acfd[:,0],2),axis=1)
        self.vectors['pwr0'].setData(pwr0)

    def applyNoise(self, level, noisemodel = np.random.randn):
        '''Apply noise of power 'noise' with model 'model' to all lag'''
        acfd = self.vectors['acfd'].getData()
        xcfd = self.vectors['xcfd'].getData()
        
        acfd += level * (noisemodel(*acfd.shape))
        xcfd += level * (noisemodel(*acfd.shape))
        
        self.vectors['acfd'].setData(acfd)
        self.vectors['xcfd'].setData(xcfd)


def main():
    test_record = rawacf_record(filename = 'sandbox/test.rawacf', scalars = DEF_SCALAR_OVERRIDES_45KM, vectors = DEF_VECTOR_OVERRIDES_45KM)

    test_record.setTime(datetime.datetime.now())
    test_record.applyNoise(.01)    
    test_record.addScatter(0, 200, 200, model = LAMBDA_FIT)

    test_record.write()
 
    test_record.close()

if __name__ == '__main__':
    main()

