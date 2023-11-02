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
from typing import Any

import mpmath as mp  # type: ignore
import numpy as np
import pytest

try:
    from rpy2 import robjects as ro  # type: ignore
    from rpy2.robjects import numpy2ri as ro_numpy2ri  # type: ignore
    from rpy2.robjects.packages import importr as ro_importr  # type: ignore
except ValueError:
    err_msg = "The R language is required to run some tests in `toms708_test.py`."
    raise RuntimeError(err_msg) from None

from pyspecials._core.toms708 import (
    ibeta,
    ibetac,
    lbeta,
    lbeta_correction,
    lgamma_difference,
)
from pyspecials._core.typing import Array

SHAPE_MAGNITUDES = ("tiny", "small", "medium", "large", "huge")

shape_cache: dict[tuple[str, np.dtype], Array] = {}


def _get_x(dtype: np.dtype) -> np.ndarray:
    epsneg = np.finfo(dtype).epsneg

    base = 10.0
    y = np.logspace(np.log10(epsneg), stop=-2, num=6, base=base, endpoint=False)
    exponent = np.floor(np.log10(y))
    tiny_values = np.ceil(y * base ** np.abs(exponent)) * base**exponent

    values = np.arange(start=0.02, stop=0.98, step=0.02)

    return np.concatenate((tiny_values, values, 1.0 - tiny_values[::-1]))


def _get_shape(magnitude: str, dtype: np.dtype, num: int = 30) -> Array:
    if magnitude in shape_cache:
        return shape_cache[(magnitude, dtype)]

    if magnitude == "tiny":
        half_num = num // 2
        epsneg = np.finfo(dtype).epsneg

        base = 10.0
        y = np.logspace(
            np.log10(epsneg),
            stop=0,
            num=num - half_num,
            base=base,
            endpoint=False,
            dtype=dtype,
        )
        exponent = np.floor(np.log10(y))
        tiny_values = np.ceil(y * base ** np.abs(exponent)) * base**exponent

        shape = np.concatenate(
            (
                tiny_values,
                np.round(
                    np.linspace(
                        start=1, stop=10, num=half_num, endpoint=False, dtype=dtype
                    ),
                    2,
                ),
            )
        )

        shape_cache[(magnitude, dtype)] = shape
        return shape

    if magnitude == "small":
        shape = np.linspace(start=10, stop=100, num=num, endpoint=False, dtype=dtype)
        shape_cache[(magnitude, dtype)] = shape
        return shape

    if magnitude == "medium":
        shape = np.linspace(start=100, stop=1_000, num=num, endpoint=False, dtype=dtype)
        shape_cache[(magnitude, dtype)] = shape
        return shape

    if magnitude == "large":
        shape = np.linspace(
            start=1_000, stop=10_000, num=num, endpoint=False, dtype=dtype
        )
        shape_cache[(magnitude, dtype)] = shape
        return shape

    if magnitude == "huge":
        shape = np.linspace(
            start=10_000, stop=100_000, num=num, endpoint=False, dtype=dtype
        )
        shape_cache[(magnitude, dtype)] = shape
        return shape

    err_msg = (
        'The argument `magnitude` must be one of `("tiny", "small", "medium", '
        '"large", "huge")`.'
    )
    raise ValueError(err_msg)


@pytest.fixture(scope="module")
def r_stats() -> Any:
    return ro_importr("stats")


@pytest.fixture(scope="module")
def numpy_converter() -> Any:
    return ro.default_converter + ro_numpy2ri.converter


@pytest.mark.parametrize(
    "a_magnitude, b_magnitude",
    itertools.product(SHAPE_MAGNITUDES, SHAPE_MAGNITUDES),
)
def test_ibeta(
    r_stats: Any, numpy_converter: Any, a_magnitude: str, b_magnitude: str
) -> None:
    space_a = _get_shape(a_magnitude, dtype=np.dtype("float64"))
    space_b = _get_shape(b_magnitude, dtype=np.dtype("float64"))
    space_x = _get_x(np.dtype("float64"))

    a, b, x = (
        np.asarray(arr)
        for arr in zip(*list(itertools.product(space_a, space_b, space_x)), strict=True)
    )

    assert a.size == (space_a.size * space_b.size * space_x.size)
    assert a.shape == b.shape
    assert a.shape == x.shape

    # Here we use the R function `stats::pbeta` because `mpmath.betainc` is very slow
    # and does not always converge. This function is based on a C translation of ACM
    # TOMS 708:
    # https://svn.r-project.org/R/trunk/src/nmath/toms708.c
    with numpy_converter.context():
        expected = r_stats.pbeta(x, a, b)

    actual = ibeta(a, b, x)

    np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=1e-20)


@pytest.mark.parametrize(
    "a_magnitude, b_magnitude",
    itertools.product(SHAPE_MAGNITUDES, SHAPE_MAGNITUDES),
)
def test_ibetac(
    r_stats: Any, numpy_converter: Any, a_magnitude: str, b_magnitude: str
) -> None:
    space_a = _get_shape(a_magnitude, dtype=np.dtype("float64"))
    space_b = _get_shape(b_magnitude, dtype=np.dtype("float64"))
    space_x = _get_x(np.dtype("float64"))

    a, b, x = (
        np.asarray(arr)
        for arr in zip(*list(itertools.product(space_a, space_b, space_x)), strict=True)
    )

    assert a.size == (space_a.size * space_b.size * space_x.size)
    assert a.shape == b.shape
    assert a.shape == x.shape

    # Here we use the R function `stats::pbeta` because `mpmath.betainc` is very slow
    # and does not always converge. This function is based on a C translation of ACM
    # TOMS 708:
    # https://svn.r-project.org/R/trunk/src/nmath/toms708.c
    with numpy_converter.context():
        expected = r_stats.pbeta(x, a, b, lower_tail=False)

    actual = ibetac(a, b, x)

    np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=1e-20)


def _mp_lbeta(a: Array, b: Array) -> Array:
    ufunc_mp_lbeta = np.frompyfunc(
        lambda x, y: mp.log(mp.beta(x, y)),
        nin=2,
        nout=1,
    )

    with mp.workdps(50):
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

    np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=0.0)


def _mp_lbeta_correction(a: Array, b: Array) -> Array:
    _a, _b = np.minimum(a, b), np.maximum(a, b)

    def _impl(x: Array, y: Array) -> Array:
        lbeta = mp.log(mp.beta(x, y))
        e = 0.5 * (mp.log(2 * mp.pi) - mp.log(y))
        h = mp.mpf(x) / mp.mpf(y)
        c = h / (1.0 + h)
        u = -(mp.mpf(x) - 0.5) * mp.log(c)
        v = mp.mpf(y) * mp.log(1.0 + h)
        return lbeta - e + u + v

    ufunc_mp_lbeta_correction = np.frompyfunc(
        lambda x, y: _impl(x, y),
        nin=2,
        nout=1,
    )

    with mp.workdps(50):
        res = ufunc_mp_lbeta_correction(_a, _b)

    return np.where((a < 8.0) | (b < 8.0), np.nan, np.asarray(res, dtype=a.dtype))


@pytest.mark.parametrize(
    "a_magnitude, b_magnitude",
    itertools.product(SHAPE_MAGNITUDES, SHAPE_MAGNITUDES),
)
def test_lbeta_correction(a_magnitude: str, b_magnitude: str) -> None:
    space_a = _get_shape(a_magnitude, dtype=np.dtype("float64"))
    space_b = _get_shape(b_magnitude, dtype=np.dtype("float64"))

    a, b = (
        np.asarray(arr)
        for arr in zip(*list(itertools.product(space_a, space_b)), strict=True)
    )

    assert a.size == (space_a.size * space_b.size)
    assert a.shape == b.shape

    expected = _mp_lbeta_correction(a, b)
    actual = lbeta_correction(a, b)

    np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=0.0)


def _mp_lgamma_difference(a: Array, b: Array) -> Array:
    ufunc_mp_lgamma_difference = np.frompyfunc(
        lambda x, y: mp.loggamma(y) - mp.loggamma(mp.mpf(x) + mp.mpf(y)),
        nin=2,
        nout=1,
    )

    with mp.workdps(50):
        res = ufunc_mp_lgamma_difference(a, b)

    return np.where(b < 8.0, np.nan, np.asarray(res, dtype=a.dtype))


@pytest.mark.parametrize(
    "a_magnitude, b_magnitude",
    itertools.product(SHAPE_MAGNITUDES, SHAPE_MAGNITUDES),
)
def test_lgamma_difference(a_magnitude: str, b_magnitude: str) -> None:
    space_a = _get_shape(a_magnitude, dtype=np.dtype("float64"))
    space_b = _get_shape(b_magnitude, dtype=np.dtype("float64"))

    a, b = (
        np.asarray(arr)
        for arr in zip(*list(itertools.product(space_a, space_b)), strict=True)
    )

    assert a.size == (space_a.size * space_b.size)
    assert a.shape == b.shape

    expected = _mp_lgamma_difference(a, b)
    actual = lgamma_difference(a, b)

    np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=0.0)
