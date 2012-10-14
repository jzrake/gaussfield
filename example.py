
import numpy as np
import matplotlib.pyplot as plt
import driving

N = 64


def example1():
    X, Y = np.mgrid[0:1:1./N, 0:1:1./N]
    field = driving.DrivingField2d(1.0, 1.0, 3, 12345)
    field.advance(1.0)
    print 'sampling field...'
    F = field.sample(X, Y)
    print 'finished'
    plt.imshow(F[:,:,0]**2 + F[:,:,1]**2)
    plt.show()


def example2():
    X, Y, Z = np.mgrid[0:1:1./N, 0:1:1./N,0:1:1./N]
    field = driving.DrivingField3d(1.0, 1.0, 3, 12345)
    field.advance(1.0)
    print 'sampling field...'
    F = field.sample(X, Y, Z)
    print 'finished'
    F2 = F[...,0]**2 + F[...,1]**2 + F[...,2]**2
    print 'mean power:', F2.mean()
    plt.imshow(F2[:,:,0])
    plt.show()


if __name__ == '__main__':
    example1()
    example2()
