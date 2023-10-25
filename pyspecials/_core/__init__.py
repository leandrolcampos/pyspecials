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

"""Implementation of the core numerical methods of PySpecials.

Please note that this module is private. All functions and objects are available in
the main `pyspecials` namespace - use that instead.

"""

from pyspecials._core.toms708 import betainc, betaincc

__all__ = [
    "betainc",
    "betaincc",
]
