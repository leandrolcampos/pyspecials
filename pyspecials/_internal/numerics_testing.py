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

"""Utilities for testing numerical accuracy."""

import numpy as np

from pyspecials._core.typing import Array, ArrayLike

__all__ = [
    "relative_error",
]


def relative_error(result: ArrayLike, truth: ArrayLike) -> Array:
    """
    Computes the relative error of `result` relative `truth`.

    The relative error is defined as `abs((result - truth) / truth)`. The computation
    of that difference and ratio are done in 64-bit precision.

    Parameters
    ----------
    result : array_like
        Array of values whose deviation to assess.
    truth : array_like
        Array of values presumed correct. Must broadcast with `result`.

    Returns
    -------
    err : array
        Float64 array of elementwise relative error values.

    """
    result = np.asarray(result, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)

    truth_is_zero = np.equal(truth, 0.0)
    safe_truth = np.where(truth_is_zero, np.ones_like(truth), truth)

    relerr = np.abs((result - truth) / safe_truth)
    relerr = np.where(truth_is_zero, np.inf, relerr)
    relerr = np.where(result == truth, np.zeros_like(truth), relerr)

    return relerr
