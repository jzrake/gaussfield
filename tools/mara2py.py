
import os
import h5py
import driving
import pickle

fname = 'SRHD-128.h5'

mara_cmd = """
json = require 'json'
h5_open_file('%s', 'r')
f = io.open('driving.dat', 'w')
f:write(json.decode(h5_read_string('driving')))
""" % fname

cmdfile = open('cmdfile.lua', 'w')
cmdfile.write(mara_cmd)
cmdfile.close()
os.system('mara cmdfile.lua')
os.system('rm cmdfile.lua')

drvdat = open('driving.dat').read()
os.system('rm driving.dat')

field = driving.DrivingField3d()
field.__setstate__(drvdat)

h5f = h5py.File(fname)

if 'new_driving' in h5f.keys():
    del h5f['new_driving']

h5f['new_driving'] = pickle.dumps(field)
