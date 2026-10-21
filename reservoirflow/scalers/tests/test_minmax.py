import unittest

import numpy as np

from reservoirflow.scalers import MinMaxScaler


class TestMinMax(unittest.TestCase):
    def test_create(self):
        MinMaxScaler((0, 1))
        MinMaxScaler(output_range=(0, 1))
        MinMaxScaler(output_range=(0, 1), input_range=None)
        MinMaxScaler(output_range=(0, 1), input_range=(10, 100))

    def test_fit(self):
        arr = np.array(
            [
                [0, 1, 1, 10, 10],
                [0, 1, 2, 20, 100],
                [0, 1, 3, 30, 1000],
                [0, 1, 4, 40, 10000],
            ]
        )
        scaler = MinMaxScaler((0, 1))
        scaler.fit(arr)
        print(scaler.transform(arr))


if __name__ == "__main__":
    unittest.main()
