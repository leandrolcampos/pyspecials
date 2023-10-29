# Copyright 2023 The PySpecials Authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for the procedures in `toms708` module."""

import itertools

import mpmath as mp  # type: ignore
import numpy as np
import pytest

from pyspecials._core.toms708 import lbeta
from pyspecials._core.typing import Array

SHAPE_MAGNITUDES = ("tiny", "small", "medium", "large", "huge")

shape_cache: dict[tuple[str, np.dtype], Array] = {}


def _get_shape(magnitude: str, dtype: np.dtype) -> Array:
    if magnitude in shape_cache:
        return shape_cache[(magnitude, dtype)]

    if magnitude == "tiny":
        epsneg = np.finfo(dtype).epsneg

        base = 10.0
        y = np.logspace(
            np.log10(epsneg), stop=0, num=20, base=base, endpoint=False, dtype=dtype
        )
        exponent = np.floor(np.log10(y))
        tiny_values = np.ceil(y * base ** np.abs(exponent)) * base**exponent

        shape = np.concatenate(
            (
                tiny_values,
                np.round(
                    np.linspace(start=1, stop=10, num=30, endpoint=False, dtype=dtype),
                    2,
                ),
            )
        )

        shape_cache[(magnitude, dtype)] = shape
        return shape

    if magnitude == "small":
        shape = np.linspace(start=10, stop=100, num=50, endpoint=False, dtype=dtype)
        shape_cache[(magnitude, dtype)] = shape
        return shape

    if magnitude == "medium":
        shape = np.linspace(start=100, stop=1_000, num=50, endpoint=False, dtype=dtype)
        shape_cache[(magnitude, dtype)] = shape
        return shape

    if magnitude == "large":
        shape = np.linspace(
            start=1_000, stop=10_000, num=50, endpoint=False, dtype=dtype
        )
        shape_cache[(magnitude, dtype)] = shape
        return shape

    if magnitude == "huge":
        shape = np.linspace(
            start=10_000, stop=100_000, num=50, endpoint=False, dtype=dtype
        )
        shape_cache[(magnitude, dtype)] = shape
        return shape

    err_msg = (
        'The argument `magnitude` must be one of `("tiny", "small", "medium", '
        '"large", "huge")`.'
    )
    raise ValueError(err_msg)


def _mp_lbeta(a: Array, b: Array) -> Array:
    ufunc_mp_lbeta = np.frompyfunc(
        lambda x, y: mp.log(mp.beta(x, y)),
        nin=2,
        nout=1,
    )

    with mp.workdps(30):
        res = ufunc_mp_lbeta(a, b)

    return np.asarray(res, dtype=a.dtype)


@pytest.mark.parametrize(
    "a_magnitude, b_magnitude",
    itertools.product(SHAPE_MAGNITUDES, SHAPE_MAGNITUDES),
)
def test_lbeta(a_magnitude: str, b_magnitude: str) -> None:
    space_a = _get_shape(a_magnitude, dtype=np.dtype("float64"))
    space_b = _get_shape(b_magnitude, dtype=np.dtype("float64"))

    a, b = (
        np.asarray(arr)
        for arr in zip(*list(itertools.product(space_a, space_b)), strict=True)
    )

    assert a.size == (space_a.size * space_b.size)
    assert a.shape == b.shape

    expected = _mp_lbeta(a, b)
    actual = lbeta(a, b)

    np.testing.assert_allclose(actual, expected, rtol=5e-14, atol=0.0)
