# script for reading in dmap strings over a socket into a dictionary
# 7/1/2015
# jon klein, jtklein@alaska.edu

import numpy as np
import pdb
import socket 
import time
import json

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
    DATACHAR:np.uint8,\
    DATASHORT:np.int16,\
    DATAINT:np.int32,\
    DATAFLOAT:np.float32,\
    DATASTRING:str}

def recv_dtype(sock, dtype, nitems = 1):
    if dtype == str:
        return recv_str(sock)
    
    dstr = ''
    while len(dstr) < dtype().nbytes * nitems:
        dstr += sock.recv(1)

    data = np.fromstring(dstr, dtype=dtype, count=nitems)
    if nitems == 1:
        return data[0]
    return data

def recv_str(sock):
    dstr = ''
    c = sock.recv(1)

    while c != NULL:
        dstr += c
        c = sock.recv(1)
    return dstr

def readPacket(sock):
    # read in header
    datacode = recv_dtype(sock, np.int32)
    # this isn't very robust.. try to find header.. or something that looks like header
    while datacode != 65537:
        datacode = recv_dtype(sock, np.int32)
    
    sze = recv_dtype(sock, np.int32)
    snum = recv_dtype(sock, np.int32)
    anum = recv_dtype(sock, np.int32)

    scalars = {}
    vectors = {}

    # read in scalars
    for s in range(snum):
        name = recv_str(sock)
        dtype = DTYPE_CODES[recv_dtype(sock, np.uint8)]
        payload = recv_dtype(sock, dtype)
        scalars[name] = payload

    # read in vectors
    for a in range(anum):
        name = recv_str(sock)
        dtype = DTYPE_CODES[recv_dtype(sock, np.uint8)]
        ndims = recv_dtype(sock, np.int32)
        dims = recv_dtype(sock, np.int32, ndims)
        payload = recv_dtype(sock, dtype, np.prod(dims))
        if ndims > 1:
            payload = np.reshape(payload, tuple(dims[::-1]))
        vectors[name] = payload
    
    return scalars, vectors

def createjson(scalars, vectors):
    json_payload = {}

    for scalar in scalars.keys():
        payload = scalars[scalar]
        if isinstance(payload, str):
            json_payload[scalar] = payload
        else:
            json_payload[scalar] = np.asscalar(payload)

    for vector in vectors.keys():
        payload = vectors[vector]
        # convert numpy array to list, recursively traverse and convert to jsonable data type..
        payload = payload.tolist()
        json_payload[vector] = payload
    return json.dumps(json_payload)


def main():
    HOST = 'superdarn.gi.alaska.edu'
    PORT = 6032

    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((HOST, PORT)) 
    while True:
        scalars, vectors = readPacket(s)
        json_str = createjson(scalars, vectors)
        print json_str

    s.close() 


if __name__ == '__main__':
    main()
