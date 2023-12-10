import unittest

import numpy as np

from reservoirflow.scalers import MinMax


class TestMinMax(unittest.TestCase):
    def test_create(self):
        MinMax((0, 1))
        MinMax(output_range=(0, 1))
        MinMax(output_range=(0, 1), input_range=None)
        MinMax(output_range=(0, 1), input_range=(10, 100))

    def test_fit(self):
        arr = np.array(
            [
                [0, 1, 1, 10, 10],
                [0, 1, 2, 20, 100],
                [0, 1, 3, 30, 1000],
                [0, 1, 4, 40, 10000],
            ]
        )
        scaler = MinMax((0, 1))
        scaler.fit(arr)
        arr_trans = scaler.transform(arr)
        arr_trans_desired = [
            [np.nan, np.nan, 0.0, 0.0, 0.0],
            [np.nan, np.nan, 0.33333333, 0.33333333, 0.00900901],
            [np.nan, np.nan, 0.66666667, 0.66666667, 0.0990991],
            [np.nan, np.nan, 1.0, 1.0, 1.0],
        ]
        np.testing.assert_almost_equal(
            arr_trans,
            arr_trans_desired,
            decimal=8,
        )


if __name__ == "__main__":
    unittest.main()
