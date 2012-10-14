
import numpy as np
import driving
import pickle
import unittest

N = 32

class TestDrivingField2d(unittest.TestCase):
    def setUp(self):
        self.seq = range(10)
        self.field = driving.DrivingField2d(1.0, 1.0, 4, 12345)
        self.X = np.zeros([N,N])
        self.Y = np.zeros([N,N])

    def test_init(self):
        F = self.field.sample(0.0, 0.0)
        self.assertEqual(F[0], 0.0)
        self.assertEqual(F[1], 0.0)

    def test_badarg(self):
        with self.assertRaises(ValueError):
            self.field.sample(self.X, 1.0)

    def test_attrib(self):
        self.assertAlmostEqual(self.field.mean_power, 1.0)
        self.assertAlmostEqual(self.field.zeta_parameter, 1.0)
        self.assertEqual(self.field.max_wavenumber, 4)

    def test_pickle(self):
        F1 = self.field.sample(0.0, 0.0)
        pk = pickle.dumps(self.field)
        field = pickle.loads(pk)
        self.assertEqual(self.field.max_wavenumber, field.max_wavenumber)
        F2 = field.sample(0.0, 0.0)
        self.assertAlmostEqual(F1[0], F2[0])
        self.assertAlmostEqual(F1[1], F2[1])


class TestDrivingField3d(unittest.TestCase):
    def setUp(self):
        self.seq = range(10)
        self.field = driving.DrivingField3d(1.0, 1.0, 4, 12345)
        self.X = np.zeros([N,N,N])
        self.Y = np.zeros([N,N,N])
        self.Z = np.zeros([N,N,N])

    def test_init(self):
        F = self.field.sample(0.0, 0.0, 0.0)
        self.assertEqual(F[0], 0.0)
        self.assertEqual(F[1], 0.0)
        self.assertEqual(F[2], 0.0)

    def test_badarg(self):
        with self.assertRaises(ValueError):
            self.field.sample(self.X, self.Y, 1.0)

    def test_attrib(self):
        self.assertAlmostEqual(self.field.mean_power, 1.0)
        self.assertAlmostEqual(self.field.zeta_parameter, 1.0)
        self.assertEqual(self.field.max_wavenumber, 4)

    def test_pickle(self):
        F1 = self.field.sample(0.0, 0.0, 0.0)
        pk = pickle.dumps(self.field)
        field = pickle.loads(pk)
        self.assertEqual(self.field.max_wavenumber, field.max_wavenumber)
        F2 = field.sample(0.0, 0.0, 0.0)
        self.assertAlmostEqual(F1[0], F2[0])
        self.assertAlmostEqual(F1[1], F2[1])
        self.assertAlmostEqual(F1[2], F2[2])


if __name__ == '__main__':
    unittest.main()
