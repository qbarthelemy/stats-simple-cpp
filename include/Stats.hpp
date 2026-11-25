#pragma once
#include <stdexcept>
#include <limits>
#include <numeric>
#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>
#include "Maths.hpp"


/**
* Statistical functions
*/
namespace Stats
{
	// --- Parametric summary statistics --- //

	/**
	* Mean.
	*
	* @tparam ContType The type of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Arithmetic mean of `x`.
	*/
	template<typename ContType>
	double mean(const ContType &x)
	{
		size_t size = x.size();
		if (size == 0)
			throw std::invalid_argument("Input has not enough values for mean.");

		double sx = std::accumulate(x.begin(), x.end(), 0.0);
		return sx / size;
	}

	/**
	* Center values with respect their mean.
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Output sequence container containing the centered version of `x`,
	* whose data type is double to keep the maximum of numerical precision.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	ContType<double, std::allocator<double>> center(const ContType<ValType, Alloc> &x)
	{
		double mean = Stats::mean(x);

		ContType<double, std::allocator<double>> x_cent(x.size());
		std::transform(x.begin(), x.end(), x_cent.begin(), std::bind2nd(std::minus<double>(), mean));
		return x_cent;
	}

	/**
	* Harmonic mean.
	*
	* @tparam ContType The type of the sequence container.
	*
	* @param x Input sequence container, containing non-zero values.
	*
	* @return Harmonic mean of `x`.
	*/
	template<typename ContType>
	double hmean(const ContType &x)
	{
		return 1.0 / Stats::mean(Maths::reciprocal(x));
	}

	/**
	* Geometric mean.
	*
	* @tparam ContType The type of the sequence container.
	*
	* @param x Input sequence container, containing strictly positive values.
	*
	* @return Geometric mean of `x`.
	*/
	template<typename ContType>
	double gmean(const ContType &x)
	{
		return std::exp(Stats::mean(Maths::log(x)));
	}

	/**
	* Power mean.
	*
	* @tparam ContType The type of the sequence container.
	*
	* @param x Input sequence container, containing strictly positive values.
	* @param exp Exponent.
	*
	* @return Power mean of `x`.
	*/
	template<typename ContType>
	double pmean(const ContType &x, double exp)
	{
		if (!Maths::is_positive(x))
			throw std::invalid_argument("Input must contain strictly positive value(s).");

		return std::pow(Stats::mean(Maths::power(x, exp)), 1.0 / exp);
	}

	/**
	* Variance.
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	* @param ddof Degree of freedom.
	*
	* @return Variance of `x`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double var(const ContType<ValType, Alloc> &x, size_t ddof = 0)
	{
		size_t size = x.size();
		if (size <= 1)
			throw std::invalid_argument("Input has not enough values for var.");
		if (size - ddof == 0)
			throw std::invalid_argument("Size minus degree of freedom is 0.");

		ContType<double, std::allocator<double>> x_cent = Stats::center(x);
		double sxx = std::inner_product(x_cent.begin(), x_cent.end(), x_cent.begin(), 0.0);
		return sxx / (size - ddof);
	}

	/**
	* Standard deviation.
	*
	* @tparam ContType The type of the sequence container.
	*
	* @param x Input sequence container.
	* @param ddof Degree of freedom.
	*
	* @return Standard deviation of `x`.
	*/
	template<typename ContType>
	double std(const ContType &x, size_t ddof = 0)
	{
		return std::sqrt(Stats::var(x, ddof));
	}

	/**
	* Harmonic standard deviation.
	*
	* @tparam ContType The type of the sequence container, containing strictly positive values.
	*
	* @param x Input sequence container.
	* @param ddof Degree of freedom.
	*
	* @return Harmonic standard deviation of `x`.
	*/
	template<typename ContType>
	double hstd(const ContType &x, size_t ddof = 0)
	{
		return 1.0 / Stats::std(Maths::reciprocal(x), ddof);
	}

	/**
	* Geometric standard deviation.
	*
	* @tparam ContType The type of the sequence container, containing strictly positive values.
	*
	* @param x Input sequence container.
	* @param ddof Degree of freedom.
	*
	* @return Geometric standard deviation of `x`.
	*/
	template<typename ContType>
	double gstd(const ContType &x, size_t ddof = 0)
	{
		return std::exp(Stats::std(Maths::log(x), ddof));
	}

	/**
	* Skewness.
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Skewness of `x`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double skew(const ContType<ValType, Alloc> &x)
	{
		size_t size = x.size();
		if (size <= 1)
			throw std::invalid_argument("Input has not enough values for skewness.");

		ContType<double, std::allocator<double>> x_cent = Stats::center(x);
		double sx2 = std::inner_product(x_cent.begin(), x_cent.end(), x_cent.begin(), 0.0);

		ContType<double, std::allocator<double>> x_cent_3 = Maths::power(x_cent, 3.0);
		double sx3 = std::accumulate(x_cent_3.begin(), x_cent_3.end(), 0.0);

		double skew = (sx3 * std::pow(size, 0.5)) / std::pow(sx2, 1.5);
		return skew;
	}

	/**
	* Kurtosis.
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Non-normalized kurtosis of `x`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double kurtosis(const ContType<ValType, Alloc> &x)
	{
		size_t size = x.size();
		if (size <= 1)
			throw std::invalid_argument("Input has not enough values for kurtosis.");

		ContType<double, std::allocator<double>> x_cent = Stats::center(x);
		double sx2 = std::inner_product(x_cent.begin(), x_cent.end(), x_cent.begin(), 0.0);

		ContType<double, std::allocator<double>> x_cent_4 = Maths::power(x_cent, 4.0);
		double sx4 = std::accumulate(x_cent_4.begin(), x_cent_4.end(), 0.0);

		double kurt = (sx4 * size) / (sx2 * sx2);
		return kurt;
	}

	// --- Nonparametric summary statistics --- //

	/**
	* Median.
	*
	* @tparam ContType The type of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Median of `x`.
	*/
	template<typename ContType>
	double median(const ContType &x)
	{
		size_t size = x.size();
		if (size == 0)
			throw std::invalid_argument("Input has not enough values for median.");

		std::vector<typename ContType::value_type> x_sort(x.begin(), x.end());
		std::sort(x_sort.begin(), x_sort.end());

		double med;
		if (size % 2 == 0)
			med = static_cast<double>(x_sort[size / 2 - 1] + x_sort[size / 2]) / 2;
		else
			med = static_cast<double>(x_sort[size / 2]);
		return med;
	}

	/**
	* Median absolute deviation.
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	* @param is_rescaled Boolean to rescale output, in order to be an estimator consistent
	* for the estimation of the standard deviation of normally distributed values.
	*
	* @return Median absolute deviation of `x`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double median_abs_deviation(const ContType<ValType, Alloc> &x, bool is_rescaled = false)
	{
		double med = Stats::median(x);

		ContType<double, std::allocator<double>> x_cent(x.size());
		std::transform(x.begin(), x.end(), x_cent.begin(), std::bind2nd(std::minus<double>(), med));
		ContType<double, std::allocator<double>> x_cent_abs = Maths::absolute(x_cent);

		double mad = Stats::median(x_cent_abs);
		if (is_rescaled)
			mad *= 1.4826;
		return mad;
	}

	// --- Transformations --- //

	/**
	* Standard score (z-score).
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	* @param ddof Degree of freedom.
	*
	* @return Output sequence container containing the z-scores of `x`,
	* whose data type is double to keep the maximum of numerical precision.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	ContType<double, std::allocator<double>> zscore(const ContType<ValType, Alloc> &x, size_t ddof = 0)
	{
		size_t size = x.size();
		if (size <= 1)
			throw std::invalid_argument("Input has not enough values for zscore.");
		if (size - ddof == 0)
			throw std::invalid_argument("Size minus degree of freedom is 0.");

		ContType<double, std::allocator<double>> x_cent = Stats::center(x);
		double sxx = std::inner_product(x_cent.begin(), x_cent.end(), x_cent.begin(), 0.0);
		double std = std::sqrt(sxx / (size - ddof));

		ContType<double, std::allocator<double>> z(size);
		std::transform(x_cent.begin(), x_cent.end(), z.begin(), std::bind2nd(std::divides<double>(), std));
		return z;
	}

	/**
	* Geometric standard score (geometric z-score).
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	* @param ddof Degree of freedom.
	*
	* @return Output sequence container containing the geometric z-scores of `x`,
	* whose data type is double to keep the maximum of numerical precision.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	ContType<double, std::allocator<double>> gzscore(const ContType<ValType, Alloc> &x, size_t ddof = 0)
	{
		return Stats::zscore(Maths::log(x), ddof);
	}

	// --- Correlation functions --- //

	/**
	* Pearson product-moment correlation coefficient.
	*
	* @tparam ContType The type of the sequence containers.
	* @tparam ValType The numeric data type of the values of the sequence containers.
	*
	* @param x, y Input sequence containers.
	*
	* @return Pearson product-moment correlation coefficient of `x` and `y`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double pearsonr(const ContType<ValType, Alloc> &x, const ContType<ValType, Alloc> &y)
	{
		size_t size = x.size();
		if (size != y.size())
			throw std::invalid_argument("Inputs have not the same size.");
		if (size == 0)
			throw std::invalid_argument("Inputs have not enough values for correlation.");

		ContType<double, std::allocator<double>> x_cent = Stats::center(x);
		ContType<double, std::allocator<double>> y_cent = Stats::center(y);

		double sxx = std::inner_product(x_cent.begin(), x_cent.end(), x_cent.begin(), 0.0);
		double syy = std::inner_product(y_cent.begin(), y_cent.end(), y_cent.begin(), 0.0);
		if (sxx <= 0 || syy <= 0)
			return std::numeric_limits<double>::quiet_NaN();

		double sxy = std::inner_product(x_cent.begin(), x_cent.end(), y_cent.begin(), 0.0);
		double r = sxy / (std::sqrt(sxx) * std::sqrt(syy));
		r = std::max(std::min(r, 1.0), -1.0);
		return r;
	}

	/**
	* Spearman rank-order correlation coefficient.
	*
	* @tparam ContType The type of the sequence containers.
	* @tparam ValType The numeric data type of the values of the sequence containers.
	*
	* @param x, y Input sequence containers.
	*
	* @return Spearman rank-order correlation coefficient of `x` and `y`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double spearmanr(const ContType<ValType, Alloc> &x, const ContType<ValType, Alloc> &y)
	{
		ContType<unsigned int, std::allocator<unsigned int>> xr = Maths::rankdata(x);
		ContType<unsigned int, std::allocator<unsigned int>> yr = Maths::rankdata(y);
		return Stats::pearsonr(xr, yr);
	}

	// --- Metrics --- //

	/**
	* Return accuracy, for binary classification.
	*
	* @tparam ContType The type of the sequence containers.
	* @tparam ValType The numeric data type of the values of the sequence containers.
	*
	* @param y_true Input sequence container containing true labels.
	* @param y_predict Input sequence container containing predicted labels.
	*
	* @return Accuracy of `y_predict` wrt. `y_true`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double accuracy_score(const ContType<ValType, Alloc> &y_true, const ContType<ValType, Alloc> &y_predict)
	{
		size_t size = y_true.size();
		if (size != y_predict.size())
			throw std::invalid_argument("Inputs have not the same size.");
		if (size == 0)
			throw std::invalid_argument("Inputs have not enough values for accuracy_score.");

		std::vector<bool, std::allocator<bool>> score(size, 0);
		std::vector<bool, std::allocator<bool>>::iterator score_it = score.begin();
		auto y_true_it = y_true.begin();
		auto y_predict_it = y_predict.begin();

		for (; y_true_it != y_true.end(); ++y_true_it, ++y_predict_it, ++score_it)
			*score_it = (static_cast<int>(*y_true_it) == static_cast<int>(*y_predict_it));
		auto c = std::count(score.begin(), score.end(), true);

		return static_cast<double>(c) / size;
	}
}