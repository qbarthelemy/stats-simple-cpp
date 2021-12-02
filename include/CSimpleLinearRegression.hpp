#include <stdexcept>
#include <limits>
#include <numeric>
#include <algorithm>
#include <vector>


/**
* Simple Linear Regression by Ordinary Least Squares.
*
* @see [sklearn.linear_model.LinearRegression](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LinearRegression.html)
*/
class SimpleLinearRegression
{
public:

	SimpleLinearRegression() { coeff = 1E42, intercept = 1E42; }

	double get_coeff() { return coeff; }

	double get_intercept() { return intercept; }

	/**
	* Fit linear model.
	*
	* @tparam ContType The type of the sequence containers.
	*
	* @param x Input sequence container containing training values.
	* @param y Input sequence container containing target values.
	*/
	template<typename ContType>
	void fit(const ContType& x, const ContType& y)
	{
		if (x.size() != y.size())
			throw std::invalid_argument("Inputs have not the same size.");
		if (x.size() < 2)
			throw std::invalid_argument("Inputs have not enough values for fit.");

		double sx = std::accumulate(x.begin(), x.end(), 0.0);
		double sy = std::accumulate(y.begin(), y.end(), 0.0);
		double sxx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
		double sxy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
		double size = static_cast<double>(x.size());
		double num = size * sxy - sx * sy;
		double denom = size * sxx - sx * sx;

		if (denom != 0)
			coeff = num / denom;
		else
			coeff = std::numeric_limits<double>::quiet_NaN();
		intercept = (sy - coeff * sx) / size;

		return;
	}

	/**
	* Predict using the linear model.
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Output sequence container containing predicted values from `x`, whose data type is double to keep the maximum of numerical precision.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	ContType<double, std::allocator<double>> predict(const ContType<ValType, Alloc>& x)
	{
		auto x_it = x.begin();
		ContType<double, std::allocator<double>> y(x.size());
		ContType<double, std::allocator<double>>::iterator y_it = y.begin();
		for (; x_it != x.end() && y_it != y.end(); ++x_it, ++y_it)
			*y_it = coeff * static_cast<double>(*x_it) + intercept;

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
	double score(const ContType<ValType, Alloc>& x, const ContType<ValType, Alloc>& y)
	{
		if (x.size() != y.size())
			throw std::invalid_argument("Inputs have not the same size.");
		if (x.size() == 0)
			throw std::invalid_argument("Inputs have not enough values for score.");

		ContType<double, std::allocator<double>> y_predict = predict(x);

		ContType<double, std::allocator<double>> res(y.size());
		std::transform(y.begin(), y.end(), y_predict.begin(), res.begin(), std::minus<double>());
		double ssres = std::inner_product(res.begin(), res.end(), res.begin(), 0.0);

		double sy = std::accumulate(y_predict.begin(), y_predict.end(), 0.0);
		double syy = std::inner_product(y_predict.begin(), y_predict.end(), y_predict.begin(), 0.0);
		double sstot = syy - sy * sy / static_cast<double>(y_predict.size());
		double r2 = 1.0 - ssres / sstot;

		return r2;
	}

protected:

	double coeff, intercept;
};
