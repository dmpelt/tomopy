# -*- coding: utf-8 -*-
import ctypes
import numpy as np

lib = ctypes.CDLL('/local/dgursoy/tomopyDev/art/art.so')
c_float_p = ctypes.POINTER(ctypes.c_float)

class cStruct(ctypes.Structure):
    _fields_ = [("numPixels", ctypes.c_int),
                ("numProjections", ctypes.c_int),
                ("numSlices", ctypes.c_int),
                ("imgN", ctypes.c_int),
                ("imgLen", ctypes.c_int),
                ("detLen", ctypes.c_int)]

class art(object):
    def __init__(self, reconInput, imgN):
        self.params = cStruct()
        self.params.numProjections = reconInput.data.shape[0]
        self.params.numSlices = reconInput.data.shape[1]
        self.params.numPixels = reconInput.data.shape[2]
        self.params.imgN = imgN;
        self.params.imgLen = reconInput.data.shape[2];
        self.params.detLen = reconInput.data.shape[2];

        # Assign/calculate values for the projection angles.
        if reconInput.angles is None:
            angles = (np.linspace(0, self.params.numProjections,
                                self.params.numProjections)
                    * 180 / self.params.numProjections).astype('float32')
        else:
            angles = np.array(reconInput.angles, dtype='float32')


        # Assign/calculate values for the centers of the slices.
        if reconInput.center is None:
            centerArr = np.ones(self.params.numSlices,
                                dtype='float32') * \
                                self.params.numPixels/2
        else:
            center = np.array(reconInput.center, dtype='float32')
            if center.size is 1:
                centerArr = np.ones(self.params.numSlices, dtype='float32') * center
            elif center.size is self.params.numSlices:
                centerArr = np.array(center, dtype='float32')
            else:
                raise ValueError('Center size must be either a scalar or equal to the number of slices.')
        print angles
        angles = np.array(angles, dtype='float32')
        center = np.array(centerArr, dtype='float32')
        rinput = np.array(reconInput.data, dtype='float32')
        self.data = np.zeros((imgN, imgN), dtype='float32')
        self.obj = lib.create(ctypes.byref(self.params),
                              angles.ctypes.data_as(c_float_p),
                              center.ctypes.data_as(c_float_p),
                              rinput.ctypes.data_as(c_float_p),
                              self.data.ctypes.data_as(c_float_p))

    def reconstruct(self):
        lib.reconstruct(self.obj)
