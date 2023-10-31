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
#
# Some of the code in this file is adapted from:
#
# numpy/numpy
# Copyright (c) 2005-2023, NumPy Developers.
# Licensed under BSD 3 clause.
#
# scipy/scipy:
# Copyright (c) 2001-2002 Enthought, Inc. 2003-2023, SciPy Developers.
# Licensed under BSD 3 clause.

"""Docstrings for ufuncs defined at `bratio.c`."""

import argparse
import os
import textwrap

MODULE_NAME = "toms708"

docdict: dict[str, str] = {}


def add_doc(name: str, doc: str) -> None:
    docdict[name] = doc


add_doc(
    "ibeta",
    r"""
    ibeta(a, b, x, out=None)

    Regularized incomplete beta function.

    Computes the regularized incomplete beta function, defined as [1]_:

    .. math::

        I_x(a, b) = \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)} \int_0^x
        t^{a-1}(1-t)^{b-1}dt,

    for :math:`0 \leq x \leq 1`.

    Parameters
    ----------
    a, b : array_like
        Positive, real-valued parameters.
    x : array_like
        Real-valued such that :math:`0 \leq x \leq 1`,
        the upper limit of integration.
    out : array, optional
        Optional output array for the function values.

    Returns
    -------
    array_like
        Value of the regularized incomplete beta function.

    See Also
    --------
    ibetac : The complement of the regularized incomplete beta function.
    lbeta :  The natural logarithm of absolute value of beta function.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/8.17

    """,
)


add_doc(
    "ibetac",
    r"""
    ibetac(a, b, x, out=None)

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
        Positive, real-valued parameters.
    x : array_like
        Real-valued such that :math:`0 \leq x \leq 1`,
        the upper limit of integration.
    out : array, optional
        Optional output array for the function values.

    Returns
    -------
    array_like
        Value of the complement the regularized incomplete beta function.

    See Also
    --------
    ibeta : The regularized incomplete beta function.
    lbeta :  The natural logarithm of absolute value of beta function.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/8.17

    """,
)


add_doc(
    "lbeta",
    r"""
    lbeta(a, b, out=None)

    Natural logarithm of beta function.

    Computes `ln(abs(beta(a, b)))` accurately.

    Parameters
    ----------
    a, b : array_like
        Positive, real-valued parameters.
    out : array, optional
        Optional output array for function values.

    Returns
    -------
    array_like
        Value of the lbeta function.

    See Also
    --------
    ibeta : The regularized incomplete beta function.
    ibetac : The complement of the regularized incomplete beta function.

    """,
)


add_doc(
    "lbeta_correction",
    r"""
    lbeta_correction(a, b, out=None)

    Error of the Stirling approximation to `lbeta(a, b)` for `a, b >= 8`.

    Parameters
    ----------
    a, b : array_like
        Real-valued parameters greater than or equal to 8.
    out : array, optional
        Optional output array for function values.

    Returns
    -------
    array_like
        Value of the lbeta_correction function.

    See Also
    --------
    lbeta :  The natural logarithm of absolute value of beta function.

    """,
)


add_doc(
    "lgamma_difference",
    r"""
    lgamma_difference(a, b, out=None)

    Difference between the natural logarithms of two gamma functions.

    Computes `ln(abs(gamma(b))) - ln(abs(gamma(a + b)))` accurately.

    Parameters
    ----------
    a, b : array_like
        Real-valued parameters.
    out : array, optional
        Optional output array for function values.

    Returns
    -------
    array_like
        Value of the lgamma_difference function.

    See Also
    --------
    lbeta :  The natural logarithm of absolute value of beta function.

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
