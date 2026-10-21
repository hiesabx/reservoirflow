import numpy as np
import pytest

from reservoirflow.scalers import MinMaxScaler


@pytest.mark.parametrize(
    "args, kwargs",
    [
        (((0, 1),), {}),
        ((), {"output_range": (0, 1)}),
        ((), {"output_range": (0, 1), "input_range": None}),
        ((), {"output_range": (0, 1), "input_range": (10, 100)}),
    ],
    ids=["positional", "output_range", "input_range_none", "input_range"],
)
def test_create(args, kwargs):
    MinMaxScaler(*args, **kwargs)


def test_fit():
    arr = np.array(
        [
            [0, 1, 1, 10, 10],
            [0, 1, 2, 20, 100],
            [0, 1, 3, 30, 1000],
            [0, 1, 4, 40, 10000],
        ]
    )
    scaler = MinMaxScaler(output_range=(0, 1))
    scaler.fit(arr, axis=0)
    arr_trans = scaler.transform(arr)
    arr_trans_desired = [
        [np.nan, np.nan, 0.0, 0.0, 0.0],
        [np.nan, np.nan, 0.33333333, 0.33333333, 0.00900901],
        [np.nan, np.nan, 0.66666667, 0.66666667, 0.0990991],
        [np.nan, np.nan, 1.0, 1.0, 1.0],
    ]
    np.testing.assert_almost_equal(arr_trans, arr_trans_desired, decimal=8)
