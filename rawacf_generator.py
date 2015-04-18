# generate synthetic rawacf files for testing fitacf against fitlomb

import numpy as np
import datetime
import sys
from pydmap_write import rawacf_record

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
DEF_SCALAR_OVERRIDES = { \
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

DEF_VECTOR_OVERRIDES = {
        'ptab' : PTAB, \
        'ltab' : LTAB, \
        'slist' : SLIST}


def main():
    dmap_file = file('temp.rawacf', 'w')
    
    test_record = rawacf_record(scalars = DEF_SCALAR_OVERRIDES, vectors = DEF_VECTOR_OVERRIDES)
    test_record.setTime(datetime.datetime.now())
    test_record.write(dmap_file)
 
    dmap_file.close()

if __name__ == '__main__':
    main()

