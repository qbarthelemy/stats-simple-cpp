# stats-simple-cpp

[![Code CppVersion](https://img.shields.io/badge/C++-11+-blue)](https://img.shields.io/badge/C++-11+-blue)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-green)](https://img.shields.io/badge/license-BSD--3--Clause-green)

Library of statistical functions and classes in simple C++,
compatible for different sequence containers of different numeric data types.

## Description

The goal of this header-only library is to mimic some functions of 
[NumPy](https://github.com/numpy/numpy) and
[SciPy](https://github.com/scipy/scipy), and some classes of
[scikit-learn](https://github.com/scikit-learn/scikit-learn) with:
- simple Python-like syntax;
- simple C++ code (compatible from C++11);
- simple dependencies (use only [C++ Standard Library](https://en.cppreference.com/w/cpp/header),
no advanced dependencies like [boost](https://github.com/boostorg/boost),
[mlpack](https://github.com/mlpack/mlpack)
or [Eigen](https://gitlab.com/libeigen/eigen));
- generic templated functions, compatible for different
[sequence containers](https://en.cppreference.com/w/cpp/named_req/SequenceContainer)
(`std::list`, `std::vector`, `std::deque`) 
of different numeric data types
(`unsigned int`, `int`, `short`, `long`, `long long`, `float`, `double`).

## Content

#### Maths.hpp

Integer input functions: `gcd`, `factorial`

Checking functions: `is_positive`

Aggregation functions: `prod`

Element-wise functions: `linear`, `absolute`, `reciprocal`,
`power`, `log`, `exp`, `sigmoid`

Others: `set`, `rankdata`

#### Stats.hpp

Summary statistics: `mean`, `center`, `hmean`, `gmean`, `pmean`,
`var`, `std`, `hstd`, `gstd`,
`skew`, `kurtosis`,
`median`, `median_abs_deviation`

Transformations: `zscore`, `gzscore`

Correlation functions: `pearsonr`, `spearmanrho`

Metrics: `accuracy_score`

#### SimpleLinearRegression.hpp

Class `SimpleLinearRegression` for linear regression,
with `fit`, `predict`, and `score` (coefficient of determination R²).

#### SimpleLogisticRegression.hpp

Class `SimpleLogisticRegression` for binary classification,
with `fit`, `predict_proba`, `predict`, and `score` (accuracy).

#### MovingZscoreDetector.hpp

Class `MovingZscoreDetector` for online anomaly detection with a moving z-score,
with `partial_fit`, and `detect`.

## Example

```
#include "SimpleLinearRegression.hpp"
#include <list>

int main()
{
	std::list<double> x = { 5.1, 13, -5, 17.5, 8 }, y = { -2, 3, 7.2, 10, -1 };
	SimpleLinearRegression slr = SimpleLinearRegression();
	slr.fit(x, y);
	x = { 42 };
	y = slr.predict(x);
}
```

## Contributing

Code must be compliant with all features listed in Description.
