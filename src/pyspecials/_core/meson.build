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

toms708_sources = ['toms708/bratio.f', 'toms708/bratio.c']

src_bratio_doc_h = custom_target('toms708_doc_generated',
  output : 'toms708_doc_generated.h',
  input : 'toms708/bratio_docstrings.py',
  command: [py, '@INPUT@', '-o', '@OUTPUT@'],
)

c_args = [
  '-DPYSPECIALS_DOC_GENERATED',
]

py.extension_module('toms708',
  [
    src_bratio_doc_h,
    toms708_sources,
  ],
  c_args: c_args,
  dependencies: py_dep,
  include_directories: np_inc,
  install: true,
  subdir: 'pyspecials/_core'
)

python_sources = [
  '__init__.py',
  'toms708.pyi',
  'typing.py',
]

py.install_sources(
  python_sources,
  subdir: 'pyspecials/_core',
  pure: false,
)
