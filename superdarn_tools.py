import numpy as np

def calc_lag_times(ltab, mplgs, mpinc):
    lags = np.float32(np.array(map(lambda x : abs(x[1]-x[0]), ltab[0:mplgs])) * (mpinc / 1e6))
    return

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

