import numpy as np


class FixedGaussianField3d(object):
    """
    Generates a Gaussian random vector field whose spatial realization is
    accomplished by a Fourier transform instead of summation of the trig
    series. It cannot be perturbed stochastically like the StochasticVectorField
    classes.
    """
    def __init__(self, size, rms=1.0, Pofk=None):
        A = np.random.uniform(-0.5, 0.5, [3] + [size]*3)
        Ak = np.zeros(A.shape, dtype=complex)
        Ks = np.zeros(A.shape, dtype=float)
        parseval_factor = (6*np.pi)**0.5 # ratio of <|Ak|^2> / <|Ax|^2>
        for a in range(3): Ak[a] = np.fft.fftn(A[a])
        Ks[0] = np.fft.fftfreq(size)[:,None,None]
        Ks[1] = np.fft.fftfreq(size)[None,:,None]
        Ks[2] = np.fft.fftfreq(size)[None,None,:]
        K2 = np.abs(Ks[0])**2 + np.abs(Ks[1])**2 + np.abs(Ks[2])**2
        K2[0,0,0] = 1.0 # prevent divide-by-zero
        Pk = Pofk(K2**0.5)
        for a in range(3): Ak[a] *= (Pk / K2)**0.5
        self._Ak = Ak
        self._Ks = Ks
        self._K2 = K2
        self._parseval_factor = parseval_factor
        self._rms = rms

    def root_mean_square(self):
        Ax = self.get_field()
        Pxs = np.abs(Ax[0])**2 + np.abs(Ax[1])**2 + np.abs(Ax[2])**2
        return Pxs.mean()**0.5

    def get_field(self, zeta=None):
        Ak = self._Ak
        Kh = self._Ks / self._K2**0.5
        if zeta:
            # ----------------------------------
            # dilatational and solenoidals parts
            # ----------------------------------
            Dk = (Ak[0]*Kh[0] + Ak[1]*Kh[1] + Ak[2]*Kh[2]) * Kh
            Sk = Ak - Dk
            Ax = np.fft.ifftn(zeta * Sk + (1.0 - zeta) * Dk).real
        else:
            Ax = np.fft.ifftn(Ak).real
        Pxs = np.abs(Ax[0])**2 + np.abs(Ax[1])**2 + np.abs(Ax[2])**2
        return self._rms * Ax / Pxs.mean()**0.5

    def power_spectrum(self, bins=128):
        Ak = self._Ak
        Ks = self._Ks
        Ps = (np.abs(Ak[0])**2 + np.abs(Ak[1])**2 + np.abs(Ak[2])**2).flatten()
        K2 = (np.abs(Ks[0])**2 + np.abs(Ks[1])**2 + np.abs(Ks[2])**2).flatten()
        vals, bins = np.histogram(K2**0.5, weights=Ps, bins=bins)
        dP = vals
        dk = 1.0 * (bins[1:] - bins[:-1])
        k0 = 0.5 * (bins[1:] + bins[:-1])
        return k0, dP/dk


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    p = -5./3.
    field = FixedGaussianField3d(128, rms=4.0, Pofk=lambda k: k**p)
    A = field.get_field()

    print field.root_mean_square()

    k, P = field.power_spectrum()

    plt.figure()
    i = np.argmax(P)
    plt.loglog(k, P)
    plt.loglog(k, P[i]*(k/k[i])**p)

    
    plt.figure()
    plt.imshow(A[0,:,:,0], interpolation='nearest')
    plt.colorbar()

    plt.show()
