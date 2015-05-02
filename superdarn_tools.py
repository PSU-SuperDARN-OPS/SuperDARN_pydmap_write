import numpy as np
import re
import subprocess
import time

C = 3e8
LAMBDA_FIT = 1
SIGMA_FIT = 2

VECTOR_SPLITTER = 'short|float|int|char'
ELEM_SPLITTER = '\n|\t'


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

    def __str__(self):
        return 'r: {}, p: {}, v: {}, w: {}, nl: {}, v_e: {}, w_e: {}'.format(self.rangegate, self.power, self.velocity, self.width, self.nlag, self.v_e, self.w_e)

def calc_lag_times(ltab, mplgs, mpinc):
    return np.float32(np.array(map(lambda x : abs(x[1]-x[0]), ltab[0:mplgs])) * (mpinc / 1e6))


def descale_velocity(v, tfreq):
    ''' 'descale' target velocity to a baseband frequency shift, v in (m/s), tfreq in kHz'''
    return v * 2. * tfreq * 1000 / C

def descale_width(w, tfreq):
    ''' 'descale' target spectral width in (m/s) to a exponential decay rate, tfreq in kHz'''

    return w * 2. * np.pi * tfreq * 1000 / C

# grabs the tail off a dmap dump of a fitacf from a remote radar 
def get_dmapdumpstring(lines = 200):
    cmd = 'dmapdump -d `ls -t ~/repos/SuperDARN_pydmap_write/sandbox/*.fitacf | head -1` | tail -n ' + str(int(lines))
    return subprocess.check_output(cmd, shell=True)

# gets data from latest dmapdump, load variables into dict, return dict 
# grab the last 200 lines or so to ensure a full record
def parse_dmapdumpstring(dumpstring):
    scandata = {}
    scan = dumpstring.split('scalars:')[-1].split('arrays:')
    scalars = scan[0].split('\n')
    vectors = re.split(VECTOR_SPLITTER, scan[1])
    for scalar in scalars:
        if scalar == '':
            continue
        assignment = scalar.split('\t')[-1].split(' = ')
        var = assignment[0].lstrip('"').rstrip('"')
        value = eval(assignment[1])
        scandata[var] = value
    for vector in vectors:
        vector = vector.split('=')

        if len(vector) <= 1:
            continue
        var = vector[0].split('"')[1]
        vecvalue = []
        for v in re.split(ELEM_SPLITTER, vector[1]):
            v = v.rstrip(',')
            if v == '':
                continue
            if v == 'inf' or v == 'nan' or v == '-nan':
                v = 'float("NaN")'
            try:
                vecvalue.append(eval(v))
            except:
                print 'error parsing vector'

        scandata[var] = np.array(vecvalue)
    return scandata

'''
def CalcNoise(pwr0, ):
    # take average of smallest ten powers at range gate 0 for lower bound on noise
    pnmin = np.mean(sorted(pwr0)[:10])

    # take 1.6 * pnmin as upper bound for noise, 
    pnmax = 1.6 * pnmin # why 1.6? because fitacf does it that way...

    noise_samples = np.array([])

    # look through good lags for ranges with pnmin, pnmax for more noise samples
    noise_ranges = (pwr0 > pnmin) * (pwr0 < pnmax)

    for r in np.nonzero(noise_ranges)[0]:
        t, samples = self._CalcSamples(r)

        noise_lags = np.nonzero((abs(samples) > pnmin) * (abs(samples) < pnmax))[0]
        noise_samples = np.append(noise_samples, abs(samples)[noise_lags])

    # set noise as average of noise samples between pnmin and pnmax
    if len(noise_samples):
        self.noise = np.mean(noise_samples)
'''

