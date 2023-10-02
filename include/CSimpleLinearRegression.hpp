#include <stdexcept>
#include <limits>
#include <numeric>
#include <algorithm>
#include <vector>


/**
* Simple linear regression, using ordinary least squares.
*
* @see [sklearn.linear_model.LinearRegression](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LinearRegression.html)
*/
class SimpleLinearRegression
{
public:

	/**
	* Create linear model.
	*/
	SimpleLinearRegression() { _coeff = 1E42, _intercept = 1E42; }

	double get_coeff() const { return _coeff; }

	double get_intercept() const { return _intercept; }

	/**
	* Fit linear model, using ordinary least squares.
	*
	* @tparam ContType The type of the sequence containers.
	*
	* @param x Input sequence container containing training values.
	* @param y Input sequence container containing target values.
	*/
	template<typename ContType>
	void fit(const ContType& x, const ContType& y)
	{
		size_t size = x.size();
		if (size != y.size())
			throw std::invalid_argument("Inputs have not the same size.");
		if (size < 2)
			throw std::invalid_argument("Inputs have not enough values for fit.");

		double sx = std::accumulate(x.begin(), x.end(), 0.0);
		double sy = std::accumulate(y.begin(), y.end(), 0.0);
		double sxx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
		double sxy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
		double num = size * sxy - sx * sy;
		double denom = size * sxx - sx * sx;

		if (denom != 0)
			_coeff = num / denom;
		else
			_coeff = std::numeric_limits<double>::quiet_NaN();
		_intercept = (sy - _coeff * sx) / size;
	}

	/**
	* Predict using the linear model.
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Output sequence container containing predicted values from `x`,
	* whose data type is double to keep the maximum of numerical precision.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	ContType<double, std::allocator<double>> predict(const ContType<ValType, Alloc>& x) const
	{
		auto x_it = x.begin();
		ContType<double, std::allocator<double>> y(x.size());
		ContType<double, std::allocator<double>>::iterator y_it = y.begin();
		for (; x_it != x.end() && y_it != y.end(); ++x_it, ++y_it)
			*y_it = _coeff * static_cast<double>(*x_it) + _intercept;

		return y;
	}

	/**
	* Return the coefficient of determination of the prediction.
	*
	* @tparam ContType The type of the sequence containers.
	* @tparam ValType The numeric data type of the values of the sequence containers.
	*
	* @param x Input sequence container containing test values.
	* @param y Input sequence container containing true predictions.
	*
	* @return R² of `predict(x)` wrt. `y`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double score(const ContType<ValType, Alloc>& x, const ContType<ValType, Alloc>& y) const
	{
		size_t size = x.size();
		if (size != y.size())
			throw std::invalid_argument("Inputs have not the same size.");
		if (size == 0)
			throw std::invalid_argument("Inputs have not enough values for score.");

		ContType<double, std::allocator<double>> y_predict = predict(x);

		ContType<double, std::allocator<double>> res(size);
		std::transform(y.begin(), y.end(), y_predict.begin(), res.begin(), std::minus<double>());
		double ssres = std::inner_product(res.begin(), res.end(), res.begin(), 0.0);

		double sy = std::accumulate(y.begin(), y.end(), 0.0);
		double syy = std::inner_product(y.begin(), y.end(), y.begin(), 0.0);
		double sstot = syy - sy * sy / size;
		double r2 = 1.0 - ssres / sstot;

		return r2;
	}

protected:

	double _coeff, _intercept;
};
