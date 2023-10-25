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

project(
  'PySpecials',
  'c', 'cython', 'fortran',
  version: '0.1.0-dev',
  license: 'Apache-2.0',
  meson_version: '>=1.2.3',
  default_options: ['warning_level=3'],
)

cc = meson.get_compiler('c')

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

if not cc.has_header('Python.h', dependencies: py_dep)
  error('Cannot compile `Python.h`. Perhaps you need to install python-dev|python-devel')
endif

np_incdir = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

np_inc = include_directories(np_incdir)

subdir('pyspecials')
