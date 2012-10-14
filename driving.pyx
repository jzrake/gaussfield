# distutils: language = c++
# distutils: sources = [src/random.cpp, src/ou-field.cpp]

import numpy as np
cimport numpy as np
from libcpp.string cimport string


cdef extern from "src/ou-field.hpp":

    cdef cppclass StochasticVectorField2d:
        void Initialize(double P0, double zeta, int k1, int seed)
        void DeSerialize(string serial)
        void AdvanceField(double dt)
        void SampleField(double x, double y, double z, double *F)
        string Serialize()
        double GetZeta()
        double GetP0()
        int GetK1()

    cdef cppclass StochasticVectorField3d:
        void Initialize(double P0, double zeta, int k1, int seed)
        void DeSerialize(string serial)
        void AdvanceField(double dt)
        void SampleField(double x, double y, double z, double *F)
        string Serialize()
        double GetZeta()
        double GetP0()
        int GetK1()


cdef class DrivingField2d(object):
    """
    Generates isotropic, 2d Gaussian random vector fields which are
    unit-correlated in time.
    """
    cdef StochasticVectorField2d *_c
    def __cinit__(self):
        self._c = new StochasticVectorField2d()

    def __dealloc__(self):
        del self._c

    def __init__(self, P0=1.0, zeta=1.0, k1=2, seed=12345):
        """
        P0: Average power in field

        zeta: Dials between 0.0 (all dilatational) and 1.0 (all solenoidal)

        k1: Largest wave-number to be included in the k-space lattice

        seed: How to seed the internal random number generator
        """
        self._c.Initialize(P0, zeta, k1, seed)

    def __reduce__(self):
        cdef string ser = self._c.Serialize()
        return (DrivingField2d, tuple(), self._c.Serialize())

    def __setstate__(self, state):
        self._c.DeSerialize(state)

    def sample(self, X, Y):
        """
        Return the field value at the coordinates X, Y. If they are arrays of
        shape [nx,ny], then return an array of shape [nx,ny,2].
        """
        cdef np.ndarray[np.double_t,ndim=1] F1
        cdef np.ndarray[np.double_t,ndim=3] F3
        cdef np.ndarray[np.double_t,ndim=2] x
        cdef np.ndarray[np.double_t,ndim=2] y
        cdef int n
        cdef double *xd, *yd
        try:
            F1 = np.zeros(2)
            self._c.SampleField(X, Y, 0.0, <double*> F1.data)
            return F1
        except: # X, Y are arrays
            try:
                assert(X.shape == Y.shape)
            except: # either AssertionError or AttributeError
                raise ValueError("coordinate arrays must be identically shaped")
            x = X
            y = Y
            xd = <double*>x.data
            yd = <double*>y.data
            F3 = np.zeros(X.shape + (2,))
            for n in range(x.size):
                self._c.SampleField(xd[n], yd[n], 0.0, <double*>F3.data + 2*n)
            return F3

    def advance(self, double dt):
        """
        Step the field forward by a time dt.
        """
        self._c.AdvanceField(dt)

    property zeta_parameter:
        def __get__(self):
            return self._c.GetZeta()
    property mean_power:
        def __get__(self):
            return self._c.GetP0()
    property max_wavenumber:
        def __get__(self):
            return self._c.GetK1()


cdef class DrivingField3d(object):
    """
    Generates isotropic, 3d Gaussian random vector fields which are
    unit-correlated in time.
    """
    cdef StochasticVectorField3d *_c
    def __cinit__(self):
        self._c = new StochasticVectorField3d()

    def __dealloc__(self):
        del self._c

    def __init__(self, P0=1.0, zeta=1.0, k1=2, seed=12345):
        """
        P0: Average power in field

        zeta: Dials between 0.0 (all dilatational) and 1.0 (all solenoidal)

        k1: Largest wave-number to be included in the k-space lattice

        seed: How to seed the internal random number generator
        """
        self._c.Initialize(P0, zeta, k1, seed)

    def __reduce__(self):
        cdef string ser = self._c.Serialize()
        return (DrivingField3d, tuple(), self._c.Serialize())

    def __setstate__(self, state):
        self._c.DeSerialize(state)

    def sample(self, X, Y, Z):
        """
        Return the field value at the coordinates X, Y, Z. If they are arrays of
        shape [nx,ny,nz], then return an array of shape [nx,ny,nz,3].
        """
        cdef np.ndarray[np.double_t,ndim=1] F1
        cdef np.ndarray[np.double_t,ndim=4] F4
        cdef np.ndarray[np.double_t,ndim=3] x
        cdef np.ndarray[np.double_t,ndim=3] y
        cdef np.ndarray[np.double_t,ndim=3] z
        cdef int n
        cdef double *xd, *yd, *zd
        try:
            F1 = np.zeros(3)
            self._c.SampleField(X, Y, Z, <double*> F1.data)
            return F1
        except: # X, Y, Z are arrays
            try:
                assert(X.shape == Y.shape and Y.shape == Z.shape)
            except: # either AssertionError or AttributeError
                raise ValueError("coordinate arrays must be identically shaped")
            x = X
            y = Y
            z = Z
            xd = <double*>x.data
            yd = <double*>y.data
            zd = <double*>z.data
            F4 = np.zeros(X.shape + (3,))
            for n in range(x.size):
                self._c.SampleField(xd[n], yd[n], zd[n], <double*>F4.data + 3*n)
            return F4

    def advance(self, double dt):
        """
        Step the field forward by a time dt.
        """
        self._c.AdvanceField(dt)

    property zeta_parameter:
        def __get__(self):
            return self._c.GetZeta()
    property mean_power:
        def __get__(self):
            return self._c.GetP0()
    property max_wavenumber:
        def __get__(self):
            return self._c.GetK1()
