# Welcome to PySpecials

PySpecials is a Python library that provides well-established numerical method implementations for selected [special functions](https://en.wikipedia.org/wiki/Special_functions).

Its primary purpose is to benchmark new numerical method implementations for these mathematical functions, using time-tested algorithms as references to assess their accuracy and execution time.

PySpecials uses the [NumPy C-API](https://numpy.org/doc/stable/reference/c-api/index.html) to seamlessly bridge the gap between Python and numerical method implementations written in Fortran, C or C++. In particular, the NumPy universal functions (UFuncs) allow PySpecials to perform element-wise operations on NumPy arrays in a highly optimized and efficient manner, supporting [vectorization](https://numpy.org/doc/stable/glossary.html#term-vectorization) and [broadcasting](https://numpy.org/doc/stable/glossary.html#term-broadcast).

## Table of Contents

- [What are special functions?](#what-are-special-functions)
- [Why the focus on special functions?](#why-the-focus-on-special-functions)
- [Installation](#installation)
- [Example Usage](#example-usage)
- [Some Implementations Available](#some-implementations-available)
  * [Beta Functions](#beta-functions)
- [Contributing](#contributing)
- [References](#references)

## What are special functions?

Special functions are particular mathematical functions that play a fundamental role in various scientific and industrial disciplines, about which many useful properties are known. They find extensive applications in physics, engineering, chemistry, computer science, and statistics, being prized for their ability to provide closed-form solutions to complex problems in these fields.

We can give some examples of special function applications in AI:

- The Gaussian Error Linear Unit (GELU) [[2](#hendrycks2016)], a high-performing neural network activation function, is defined based on the [Gauss error](https://en.wikipedia.org/wiki/Error_function) function.

- Using numerical methods for [Bessel](https://en.wikipedia.org/wiki/Bessel_function), [incomplete beta](https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function), and [incomplete gamma](https://en.wikipedia.org/wiki/Incomplete_gamma_function) functions, we can implicitly differentiate [[1](#figurnov2018)] cumulative distribution functions that are expressed in terms of these special functions, and then train probabilistic models with, for instance, von Mises, gamma, and beta latent variables.

## Why the focus on special functions?

Beyond the practical importance of special functions in scientific and industrial applications, finding accurate and efficient ways to work with them can be an enjoyable brain-teaser for those who love math and computer science.

## Installation

Currently, the only way to install PySpecials is building it from source with [Meson](https://mesonbuild.com/). This requires Python 3.11 or newer.

**We strongly recommend using a virtual environment to isolate PySpecials' dependencies from your system's environment.**

First, clone the repository and install the building tools:

```
git clone https://github.com/leandrolcampos/pyspecials.git
cd pyspecials
python -m pip install -r build_requirements.txt
```

Then install PySpecials:

```
python -m pip install .
```

If you want to test PySpecials, you also need to install the R language in your system. PySpecials is tested with R version 4.3.1. After cloning the repository and installing the building tools and the R language, execute the following commands to install PySpecials and test it:

```
python -m pip install --editable ".[test]"
pytest
```

## Example Usage

The following code snippet shows how to compute the regularized incomplete beta function for given broadcastable NumPy arrays:

```python
>>> import numpy as np
>>> import pyspecials as ps
>>> a = np.array([0.10, 1.00, 10.0])
>>> b = np.array([[0.30], [3.00], [30.0]])
>>> x = np.array([0.01, 0.50, 0.99])
>>> ps.ibeta(a, b, x)
array([[0.49207421, 0.1877476 , 0.45864671],
       [0.72743419, 0.875     , 0.99979438],
       [0.90766467, 1.        , 1.        ]])
```

## Some Implementations Available

### Beta Functions

| Function | Description | Reference |
|----------|-------------|-----------|
| `ibeta(a, b, x[, out])` | Regularized incomplete beta function | [ACM TOMS 708](https://dl.acm.org/doi/10.1145/131766.131776) |
| `ibetac(a, b, x[, out])` | Complement of the regularized incomplete beta function | [ACM TOMS 708](https://dl.acm.org/doi/10.1145/131766.131776) |
| `lbeta(a, b[, out])` | Natural logarithm of absolute value of beta function | [ACM TOMS 708](https://dl.acm.org/doi/10.1145/131766.131776) |

## Contributing

We are not accepting pull requests at this time. However, you can contribute by reporting issues or suggesting features through the creation of a GitHub issue [here](https://github.com/leandrolcampos/pyspecials/issues).

## References

<a id="figurnov2018">[1]</a>
Figurnov, Mikhail, Shakir Mohamed, and Andriy Mnih. "Implicit reparameterization gradients." _Advances in neural information processing systems_ 31 (2018). [[Link](https://arxiv.org/abs/1805.08498)]

<a id="hendrycks2016">[2]</a>
Hendrycks, Dan, and Kevin Gimpel. "Gaussian error linear units (gelus)." _arXiv preprint arXiv:1606.08415_ (2016). [[Link](https://arxiv.org/abs/1606.08415)]