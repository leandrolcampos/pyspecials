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

"""Tests for the numerics testing utilities."""

import numpy as np

from pyspecials._internal.numerics_testing import relative_error


def test_relative_error() -> None:
    result = np.array([0.0, 1.1, 2.0, 3.0])
    truth = np.array([0.0, 1.0, 0.0, 3.0])
    expected = np.array([0.0, 0.1, np.inf, 0.0])

    actual = relative_error(result, truth)
    np.testing.assert_allclose(actual, expected, rtol=1e-15, atol=0.0)
