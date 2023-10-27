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

"""Docstrings for ufuncs defined at `bratio.c`."""

import argparse
import os
import textwrap

MODULE_NAME = "tosm708"

docdict: dict[str, str] = {}


def add_doc(name: str, doc: str) -> None:
    docdict[name] = doc


add_doc(
    "betainc",
    r"""
    betainc(a, b, x, out=None)

    Regularized incomplete beta function.

    Computes the regularized incomplete beta function, defined as [1]_:

    .. math::

        I_x(a, b) = \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)} \int_0^x
        t^{a-1}(1-t)^{b-1}dt,

    for :math:`0 \leq x \leq 1`.

    This function is the cumulative distribution function for the beta
    distribution; its range is [0, 1].

    Parameters
    ----------
    a, b : array_like
           Positive, real-valued parameters
    x : array_like
        Real-valued such that :math:`0 \leq x \leq 1`,
        the upper limit of integration
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Value of the regularized incomplete beta function

    See Also
    --------
    betaincc : complement of the regularized incomplete beta function

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/8.17

    """,
)


add_doc(
    "betaincc",
    r"""
    betaincc(a, b, x, out=None)

    Complement of the regularized incomplete beta function.

    Computes the complement of the regularized incomplete beta function,
    defined as [1]_:

    .. math::

        \bar{I}_x(a, b) = 1 - I_x(a, b)
                        = 1 - \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)} \int_0^x
                                  t^{a-1}(1-t)^{b-1}dt,

    for :math:`0 \leq x \leq 1`.

    Parameters
    ----------
    a, b : array_like
           Positive, real-valued parameters
    x : array_like
        Real-valued such that :math:`0 \leq x \leq 1`,
        the upper limit of integration
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Value of the regularized incomplete beta function

    See Also
    --------
    betainc : regularized incomplete beta function

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/8.17

    """,
)


def normalize_doc(doc: str) -> str:
    doc = textwrap.dedent(doc).strip()
    doc = doc.encode("unicode-escape").decode("ascii")
    doc = doc.replace(r'"', r"\"")
    doc = doc.replace(r"'", r"\'")
    # Split the docstring because some compilers (like MS) do not like big
    # string literal in C code. We split at endlines because textwrap.wrap
    # do not play well with \n
    doc = '\\n""'.join(doc.split(r"\n"))
    return doc


def write_code(outfile: str) -> None:
    flag = f"PYSPECIALS__CORE_INCLUDE_{MODULE_NAME.upper()}_DOC_GENERATED_H_"
    with open(outfile, "w") as fid:
        fid.write(f"#ifndef {flag}\n#define {flag}\n")
        for name, doc in docdict.items():
            cdef_name = f"DOC_{name}"
            cdef_doc = normalize_doc(doc)
            fid.write(f'#define {cdef_name} "{cdef_doc}"\n')
        fid.write(f"#endif //{flag}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o", "--outfile", type=str, help="Relative path to the output file"
    )
    args = parser.parse_args()

    outfile = os.path.join(os.getcwd(), args.outfile)
    write_code(outfile)


if __name__ == "__main__":
    main()
