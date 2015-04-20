# generate synthetic rawacf files for testing fitacf against fitlomb

import numpy as np
import datetime
import sys
from pydmap_write import dmap_record
from superdarn_tools import calc_lag_times 

NAVE = 30
LAGFR = 1200
SMSEP = 300
NOISE_SEARCH = 5
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
    def __init__(self, scalars = {}, vectors = {}):
        dmap_record.__init__(self)

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
        self.addVectorBlank('acfd', 0, (self.scalars['nrang'].getData(), self.scalars['mplgs'].getData(), 2), np.int16)
        self.addVectorBlank('xcfd', 0, (self.scalars['nrang'].getData(), self.scalars['mplgs'].getData(), 2), np.int16)

        for v in vectors:
            self.vectors[v].setData(vectors[v])
        
    def setDefaults(rsep = 45):
        ''' sets sane default parameters for a given rsep (km)'''
        if rsep == 45:
            self.setData(scalars = DEF_SCALAR_OVERRIDES_45KM, vectors = DEF_VECTOR_OVERRIDES_45KM)
        else:
            NotImplementedError('rsep {} does not have defaults'.format(rsep)) 
    
    def addScatter(self, rgate, velocity, spectral_width):
        ''' add scatter to acfd at gate rgate with velocity (m/s) and spectral_width (m/s) '''
        lagt = superdarn_tools.calc_lag_times(self.vectors['ltab'].data, self.scalars['mplgs'].data, self.scalars['mpinc'].data)
        # convert velocity to phase change per mpinc?
        # convert spectral width to something useful something
        #envelope = np.exp(-np.pow(a

        '''
                    double a = alfs[alfidx];
                    double env = exp(-pow(a, (double) env_model) * powf(times_secs[t], env_model)) * lagmask[t]; 
                    float vm = vlos_mps[vidx] * sqrtf(1 -(powf(fcrit_hz[fcritidx],2)/powf(tfreq_hz[t],2)));
                    float theta = 2 * PI * times_secs[t] * 2 * tfreq_hz[t] * (vm / C_MPS);
                    
                    #if USESHARED
                        r += env * (cosf(theta) * s_isamples[t] + sinf(theta) * s_qsamples[t]);
                        i += env * (sinf(theta) * s_isamples[t] - cosf(theta) * s_qsamples[t]);
                    #else
                        r += env * (cosf(theta) * isamples[t] + sinf(theta) * qsamples[t]);
                        i += env * (sinf(theta) * isamples[t] - cosf(theta) * qsamples[t]);
                    #endif

                    envsum += pow(env, 2);
                }
        '''
    def addTarget(self, rgate, velocity, spectral_width):
        ''' add targets to a list of targets to generate scatter from using generateScatter '''
        pass

    def generateScatter(self):
        ''' generate scatter including cross-range interference using a list of targets '''
        pass

    def applyNoise(self, level, noisemodel = np.random.randn):
        '''Apply noise of power 'noise' with model 'model' to all lag'''
        acfd = self.vectors['acfd'].getData()
        xcfd = self.vectors['xcfd'].getData()
        
        acfd += level * (noisemodel(*acfd.shape) + 1j * noisemodel(*acfd.shape))
        xcfd += level * (noisemodel(*acfd.shape) + 1j * noisemodel(*acfd.shape))
                
        self.vectors['acfd'].setData(acfd)
        self.vectors['xcfd'].setData(xcfd)


        

def main():
    dmap_file = file('temp.rawacf', 'w')
    test_record = rawacf_record(scalars = DEF_SCALAR_OVERRIDES_45KM, vectors = DEF_VECTOR_OVERRIDES_45KM)

    test_record.setTime(datetime.datetime.now())

    test_record.write(dmap_file)
 
    dmap_file.close()

if __name__ == '__main__':
    main()

