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

"""Static type annotations."""

import numpy as np

__all__ = [
    "Array",
    "ArrayLike",
]


# Annotation for NumPy array.
Array = np.ndarray

# Annotation for any value that is safe to implicitly cast to a NumPy array.
ArrayLike = (
    Array
    | bool
    | int
    | float
    | np.bool_
    | np.int32
    | np.int64
    | np.float32
    | np.float64
)
