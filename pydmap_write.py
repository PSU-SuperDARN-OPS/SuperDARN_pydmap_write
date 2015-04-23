# script for generating dmap files in python
# jon klein, jtklein@alaska.edu

import numpy as np
import collections
import struct
import pdb

DATACODE = 33
DATACHAR = 1
DATASHORT = 2
DATAINT = 3
DATAFLOAT = 4
DATADOUBLE = 8
DATASTRING = 9
DATALONG = 10
DATAUCHAR = 16
DATAUSHORT = 17
DATAUINT = 18
DATAULONG = 19
DATAMAP = 255
NULL = chr(0)

DTYPE_CODES = { \
    np.uint8:DATACHAR,\
    np.int16:DATASHORT,\
    np.int32:DATAINT,\
    np.float32:DATAFLOAT,\
    str:DATASTRING}

class dmap_var:
    def __init__(self, data, dtype):
        self.setType(dtype)
        self.setData(data)

    def setData(self, data):
        if type(data) == np.ndarray:
            self.data = np.array(data, dtype = self.dtype)
        else:
            self.data = self.dtype(data)

    def setType(self, dtype):
        self.dtype = dtype

    def getData(self):
        return self.data

    def getDmapPack(self, name):
        datastr = name + NULL + chr(DTYPE_CODES[self.dtype])
        if type(self.data) == np.ndarray:
            dims = self.data.shape[::-1]

            datastr += np.int32(len(dims)).tobytes()
            for d in dims:
                datastr += np.int32(d).tobytes()

        if self.dtype != str:
            datastr += self.data.tobytes()
        else:
            datastr += (self.data + NULL)

        return datastr

class dmap_record(object):
    def __init__(self, filename = ''):
        self.scalars = collections.OrderedDict()
        self.vectors = collections.OrderedDict()
        self.filename = filename

        self.addScalar('radar.revision.major', 0, np.uint8)
        self.addScalar('radar.revision.minor', 0, np.uint8)
        self.addScalar('origin.code', 0, np.uint8)
        self.addScalar('origin.time', 'test time string', str)
        self.addScalar('origin.command', 'pydmap_write.py', str)

        self.addScalar('cp', 0, np.int16)
        self.addScalar('stid', 0, np.int16)
        self.addScalar('time.yr', 0, np.int16)
        self.addScalar('time.mo', 0, np.int16)
        self.addScalar('time.dy', 0, np.int16)
        self.addScalar('time.hr', 0, np.int16)
        self.addScalar('time.mt', 0, np.int16)
        self.addScalar('time.sc', 0, np.int16)
        self.addScalar('time.us', 0, np.int32)
        self.addScalar('txpow', 0, np.int16)
        self.addScalar('nave', 0, np.int16)
        self.addScalar('atten', 0, np.int16)
        self.addScalar('lagfr', 0, np.int16)
        self.addScalar('smsep', 0, np.int16)
        self.addScalar('ercod', 0, np.int16)
        self.addScalar('stat.agc', 0, np.int16)
        self.addScalar('noise.search', 0, np.float32)
        self.addScalar('noise.mean', 0, np.float32)
        self.addScalar('channel', 0, np.int16)
        self.addScalar('bmnum', 0, np.int16)
        self.addScalar('bmazm', 0, np.float32)
        self.addScalar('scan', 0, np.int16)
        self.addScalar('offset', 0, np.int16)
        self.addScalar('rxrise', 0, np.int16)
        self.addScalar('intt.sc', 0, np.int16)
        self.addScalar('intt.us', 0, np.int32)
        self.addScalar('txpl', 0, np.int16)
        self.addScalar('mpinc', 0, np.int16)
        self.addScalar('mppul', 8, np.int16)
        self.addScalar('mplgs', 3, np.int16)
        self.addScalar('nrang', 4, np.int16)
        self.addScalar('frang', 0, np.int16)
        self.addScalar('rsep', 0, np.int16)
        self.addScalar('xcf', 0, np.int16)
        self.addScalar('tfreq', 0, np.int16)
        self.addScalar('mxpwr', 0, np.int32)
        self.addScalar('lvmax', 0, np.int32)
        self.addScalar('combf', 0, str)
        
        if self.filename:
            self.dmap_file = file(self.filename, 'w')

    def setData(self, scalars, vectors):
        for s in scalars:
            self.scalars[s].setData(scalars[s])
        for v in vectors:
            self.scalars[s].setData(scalars[s])

    # sets time from input datetime dt 
    def setTime(self, dt):
        self.scalars['time.yr'].setData(dt.year)
        self.scalars['time.mo'].setData(dt.month)
        self.scalars['time.dy'].setData(dt.day)
        self.scalars['time.hr'].setData(dt.hour)
        self.scalars['time.mt'].setData(dt.minute)
        self.scalars['time.sc'].setData(dt.second)
        self.scalars['time.us'].setData(dt.microsecond)


    def _makeHeader(self):
        datacode = DATACODE
        sze = len(self.vector_str) + len(self.scalar_str) + 4 * np.int32().nbytes
        anum = len(self.vectors)
        snum = len(self.scalars)

        self.header = struct.pack('iiii', datacode, sze, snum, anum)
    
    def _makeScalarStr(self):
        self.scalar_str = ''
        for s in self.scalars:
            self.scalar_str += self.scalars[s].getDmapPack(s)
   
    def _makeVectorStr(self):
        self.vector_str = ''
        for v in self.vectors:
            self.vector_str += self.vectors[v].getDmapPack(v)
   
    def write(self):
        self._makeScalarStr()
        self._makeVectorStr()
        self._makeHeader()
        
        self.dmap_file.write(self.header)
        self.dmap_file.write(self.scalar_str)
        self.dmap_file.write(self.vector_str)

    def close(self):
        if self.filename:
            self.dmap_file.close()

    def addScalar(self, name, val, dtype):
        self.scalars[name] = dmap_var(val, dtype)
    
    # create an array of type dtype of shape shape with default value dval
    def addVectorBlank(self, name, dval, shape, dtype):
        val = dval * np.ones(shape, dtype=dtype) 
        self.vectors[name] = dmap_var(val, dtype)

    def addVector(self, name, data, dtype):
        self.addVectorBlank(name, 0, data.shape, dtype)
        self.vectors['name'].setData(data)


class fitacf_record(dmap_record):
    def __init__(self, scalars = {}, vectors = {}):
        dmap_record.__init__(self) 
        self.addScalar('fitacf.revision.major', 0, np.int32)
        self.addScalar('fitacf.revision.minor', 0, np.int32)
        self.addScalar('noise.sky', 0, np.float32)
        self.addScalar('noise.lag0', 0, np.float32)
        self.addScalar('noise.vel', 0, np.float32)

        for s in scalars:
            self.scalars[s].setData(scalars[s])

        self.addVectorBlank('ptab', 0, (self.scalars['mppul'].getData()), np.int16)
        self.addVectorBlank('ltab', 0, (2, self.scalars['mplgs'].getData()), np.int16)
        self.addVectorBlank('pwr0', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('slist', 0, (self.scalars['nrang'].getData()), np.int16)
        self.addVectorBlank('nlag', 0, (self.scalars['nrang'].getData()), np.int16)
        self.addVectorBlank('qflg', 0, (self.scalars['nrang'].getData()), np.uint8)
        self.addVectorBlank('gflg', 0, (self.scalars['nrang'].getData()), np.uint8)
        self.addVectorBlank('p_l', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('p_l_e', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('p_s', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('p_s_e', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('v', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('v_e', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('w_l', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('w_l_e', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('w_s', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('w_s_e', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('sd_l', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('sd_s', 0, (self.scalars['nrang'].getData()), np.float32)
        self.addVectorBlank('sd_phi', 0, (self.scalars['nrang'].getData()), np.float32)
        
        for v in vectors:
            self.vectors[v].setData(vectors[v])


def main():
    dmap_file = file('sandbox/temp.rawacf', 'w')
    dmap_r = dmap_record()
   
    dmap_r.write(dmap_file)

    dmap_file.close()


if __name__ == '__main__':
    main()
