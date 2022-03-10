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
	double mean(const ContType& x)
	{
		size_t size = x.size();
		if (size == 0)
			throw std::invalid_argument("Input has not enough values for mean.");

		double sx = std::accumulate(x.begin(), x.end(), 0.0);
		return sx / size;
	}

	/**
	* Harmonic mean.
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container, containing non-zero values.
	*
	* @return Harmonic mean of `x`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double hmean(const ContType<ValType, Alloc>& x)
	{
		size_t size = x.size();
		if (size == 0)
			throw std::invalid_argument("Input has not enough values for hmean.");

		ContType<double, std::allocator<double>> x_inv = Maths::reciprocal(x);
		double sxinv = std::accumulate(x_inv.begin(), x_inv.end(), 0.0);
		return size / sxinv;
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
	double gmean(const ContType& x)
	{
		size_t size = x.size();
		if (size == 0)
			throw std::invalid_argument("Input has not enough values for gmean.");
		if (!Maths::is_positive(x))
			throw std::invalid_argument("Input contains negative value(s).");

		return std::pow(Maths::prod(x), 1.0 / size);
	}

	/**
	* Power mean.
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container, containing strictly positive values.
	* @param exp Exponent.
	*
	* @return Power mean of `x`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double pmean(const ContType<ValType, Alloc>& x, double exp)
	{
		size_t size = x.size();
		if (size == 0)
			throw std::invalid_argument("Input has not enough values for pmean.");
		if (!Maths::is_positive(x))
			throw std::invalid_argument("Input contains negative value(s).");

		ContType<double, std::allocator<double>> x_pow = Maths::power(x, exp);
		double sxpow = std::accumulate(x_pow.begin(), x_pow.end(), 0.0);
		return std::pow(sxpow / size, 1.0 / exp);
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
	double var(const ContType<ValType, Alloc>& x, size_t ddof = 0)
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
	double std(const ContType& x, size_t ddof = 0)
	{
		return std::sqrt(Stats::var(x, ddof));
	}

	/**
	* Harmonic standard deviation.
	*
	* @tparam ContType The type of the sequence container, containing strictly positive values.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	* @param ddof Degree of freedom.
	*
	* @return Harmonic standard deviation of `x`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double hstd(const ContType<ValType, Alloc>& x, size_t ddof = 0)
	{
		ContType<double, std::allocator<double>> x_inv = Maths::reciprocal(x);
		return 1.0 / Stats::std(x_inv, ddof);
	}

	/**
	* Geometric standard deviation.
	*
	* @tparam ContType The type of the sequence container, containing strictly positive values.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	* @param ddof Degree of freedom.
	*
	* @return Geometric standard deviation of `x`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double gstd(const ContType<ValType, Alloc>& x, size_t ddof = 0)
	{
		ContType<double, std::allocator<double>> x_log = Maths::log(x);
		return std::exp(Stats::std(x_log, ddof));
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
	double skewness(const ContType<ValType, Alloc>& x)
	{
		size_t size = x.size();
		if (size <= 1)
			throw std::invalid_argument("Input has not enough values for skewness.");

		ContType<double, std::allocator<double>> x_cent = Stats::center(x);
		double sx2 = std::inner_product(x_cent.begin(), x_cent.end(), x_cent.begin(), 0.0);

		ContType<double, std::allocator<double>> x_cent_3 = Maths::pow(x_cent, 3.0);
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
	double kurtosis(const ContType<ValType, Alloc>& x)
	{
		size_t size = x.size();
		if (size <= 1)
			throw std::invalid_argument("Input has not enough values for kurtosis.");

		ContType<double, std::allocator<double>> x_cent = Stats::center(x);
		double sx2 = std::inner_product(x_cent.begin(), x_cent.end(), x_cent.begin(), 0.0);

		ContType<double, std::allocator<double>> x_cent_4 = Maths::pow(x_cent, 4.0);
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
	double median(const ContType& x)
	{
		size_t size = x.size();
		if (size == 0)
			throw std::invalid_argument("Input has not enough values for median.");

		std::vector<ContType::value_type> x_sort(x.begin(), x.end());
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
	double median_abs_deviation(const ContType<ValType, Alloc>& x, bool is_rescaled = false)
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
	typename ContType<double, std::allocator<double>> center(const ContType<ValType, Alloc>& x)
	{
		double mean = Stats::mean(x);

		ContType<double, std::allocator<double>> x_cent(x.size());
		std::transform(x.begin(), x.end(), x_cent.begin(), std::bind2nd(std::minus<double>(), mean));
		return x_cent;
	}

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
	typename ContType<double, std::allocator<double>> zscore(const ContType<ValType, Alloc>& x, size_t ddof = 0)
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
	typename ContType<double, std::allocator<double>> gzscore(const ContType<ValType, Alloc>& x, size_t ddof = 0)
	{
		ContType<double, std::allocator<double>> x_log = Maths::log(x);
		ContType<double, std::allocator<double>> gz = Stats::zscore(x_log, ddof);
		return gz;
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
	double pearsonr(const ContType<ValType, Alloc>& x, const ContType<ValType, Alloc>& y)
	{
		size_t size = x.size();
		if (size != y.size())
			throw std::invalid_argument("Inputs have not the same size.");
		if (size == 0)
			throw std::invalid_argument("Inputs have not enough values for pearsonr.");

		ContType<double, std::allocator<double>> x_cent = Stats::center(x);
		ContType<double, std::allocator<double>> y_cent = Stats::center(y);

		double sxx = std::inner_product(x_cent.begin(), x_cent.end(), x_cent.begin(), 0.0);
		double syy = std::inner_product(y_cent.begin(), y_cent.end(), y_cent.begin(), 0.0);
		double sxy = std::inner_product(x_cent.begin(), x_cent.end(), y_cent.begin(), 0.0);

		double denom = std::sqrt(sxx) * std::sqrt(syy);
		if (denom <= 0)
			return std::numeric_limits<double>::quiet_NaN();

		double r = sxy / denom;
		r = max(min(r, 1.0), -1.0);
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
	double spearmanr(const ContType<ValType, Alloc>& x, const ContType<ValType, Alloc>& y)
	{
		size_t size = x.size();
		if (size != y.size())
			throw std::invalid_argument("Inputs have not the same size.");
		if (size == 0)
			throw std::invalid_argument("Inputs have not enough values for spearmanr.");

		ContType<unsigned int, std::allocator<unsigned int>> xr = Stats::rankdata(x);
		ContType<unsigned int, std::allocator<unsigned int>> yr = Stats::rankdata(y);

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
	double accuracy_score(const ContType<ValType, Alloc>& y_true, const ContType<ValType, Alloc>& y_predict)
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

	// --- Others --- //

	template<typename ValType>
	bool op_sort_increase(const std::pair<ValType, unsigned int>& x1, const std::pair<ValType, unsigned int>& x2)
	{
		return x1.first < x2.first;
	}

	template<typename ValType>
	bool op_sort_decrease(const std::pair<ValType, unsigned int>& x1, const std::pair<ValType, unsigned int>& x2)
	{
		return x1.first > x2.first;
	}

	/**
	* Assign ranks to data.
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Output sequence container containing increasing ranks of values of `x`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	typename ContType<unsigned int, std::allocator<unsigned int>> rankdata(const ContType<ValType, Alloc>& x)
	{
		size_t size = x.size(), i = 0;

		auto x_it = x.begin();
		std::vector<std::pair<ValType, unsigned int>> xr(size);
		std::vector<std::pair<ValType, unsigned int>>::iterator xr_it = xr.begin();
		for (i = 0; x_it != x.end() && xr_it != xr.end(); ++i, ++x_it, ++xr_it)
			*xr_it = std::make_pair(*x_it, static_cast<unsigned int>(i));

		std::sort(xr.begin(), xr.end(), Stats::op_sort_increase<ValType>);

		ContType<unsigned int, std::allocator<unsigned int>> ranks(size);
		ContType<unsigned int, std::allocator<unsigned int>>::iterator r_it = ranks.begin();
		for (xr_it = xr.begin(); xr_it != xr.end() && r_it != ranks.end(); ++xr_it, ++r_it)
			*r_it = (*xr_it).second;

		return ranks;
	}

}