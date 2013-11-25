# -*- coding: utf-8 -*-
import ctypes
import numpy as np

lib = ctypes.CDLL('/local/dgursoy/tomopyDev/art/art.so')
c_float_p = ctypes.POINTER(ctypes.c_float)

class cStruct(ctypes.Structure):
    _fields_ = [("numPixels", ctypes.c_int),
                ("numProjections", ctypes.c_int),
                ("numSlices", ctypes.c_int),
                ("imgN", ctypes.c_int)]

class art(object):
    def __init__(self, reconInput, imgN):
        self.params = cStruct()
        self.params.numProjections = reconInput.data.shape[0]
        self.params.numSlices = reconInput.data.shape[0]
        self.params.numPixels = reconInput.data.shape[0]
        self.params.imgN = imgN;

        # Assign/calculate values for the projection angles.
        if reconInput.angles is None:
            angles = (np.linspace(0, self.params.numProjections,
                                self.params.numProjections)
                    * 180 / self.params.numProjections).astype('float32')

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
        print " "
        print center
        angles = np.array(angles, dtype='float32')
        center = np.array(centerArr, dtype='float32')
        rinput = np.ones((4, 6, 5), dtype='float32')
        self.routput = np.empty((6, 5, 5), dtype='float32')
        self.obj = lib.create(ctypes.byref(self.params),
                              angles.ctypes.data_as(c_float_p),
                              center.ctypes.data_as(c_float_p),
                              rinput.ctypes.data_as(c_float_p),
                              self.routput.ctypes.data_as(c_float_p))

    def reconstruct(self):
        lib.reconstruct(self.obj)
