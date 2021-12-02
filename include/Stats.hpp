#include <stdexcept>
#include <limits>
#include <numeric>
#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>


namespace Stats
{
	// --- SummaryStatistics --- //

	/**
	* Mean.
	*
	* @tparam ContType The type of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Mean of `x`.
	*/
	template<typename ContType>
	double mean(const ContType& x)
	{
		if (x.size() == 0)
			throw std::invalid_argument("Input has not enough values for mean.");

		double sum = std::accumulate(x.begin(), x.end(), 0.0);
		double mean = sum / static_cast<double>(x.size());
		return mean;
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
		if (x.size() <= 1)
			throw std::invalid_argument("Input has not enough values for var.");
		if (x.size() - ddof == 0)
			throw std::invalid_argument("Size minus degree of freedom is 0.");

		ContType<double, std::allocator<double>> x_centered = Stats::center(x);
		double sxx = std::inner_product(x_centered.begin(), x_centered.end(), x_centered.begin(), 0.0);
		double var = sxx / static_cast<double>(x.size() - ddof);
		return var;
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
		double std = std::sqrt(Stats::var(x, ddof));
		return std;
	}

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

		std::vector<ContType::value_type> x_sorted(x.begin(), x.end());
		std::sort(x_sorted.begin(), x_sorted.end());

		double med;
		if (size % 2 == 0)
			med = static_cast<double>(x_sorted[size / 2 - 1] + x_sorted[size / 2]) / 2;
		else
			med = static_cast<double>(x_sorted[size / 2]);
		return med;
	}

	/**
	* Median absolute deviation.
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Median absolute deviation of `x`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double median_abs_deviation(const ContType<ValType, Alloc>& x)
	{
		double med = Stats::median(x);

		ContType<double, std::allocator<double>> x_centered(x.size()), x_centered_abs(x.size());
		std::transform(x.begin(), x.end(), x_centered.begin(), std::bind2nd(std::minus<double>(), med));
		std::transform(x_centered.begin(), x_centered.end(), x_centered_abs.begin(), static_cast<double(*)(double)>(&std::abs));

		double mad = Stats::median(x_centered_abs);
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
	* @return Output sequence container containing the centered version of `x`, whose data type is double to keep the maximum of numerical precision.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	typename ContType<double, std::allocator<double>> center(const ContType<ValType, Alloc>& x)
	{
		double mean = Stats::mean(x);

		ContType<double, std::allocator<double>> x_centered(x.size());
		std::transform(x.begin(), x.end(), x_centered.begin(), std::bind2nd(std::minus<double>(), mean));
		return x_centered;
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
	* @return Output sequence container containing the z-scores of `x`, whose data type is double to keep the maximum of numerical precision.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	typename ContType<double, std::allocator<double>> zscore(const ContType<ValType, Alloc>& x, size_t ddof = 0)
	{
		size_t size = x.size();
		if (size <= 1)
			throw std::invalid_argument("Input has not enough values for zscore.");
		if (size - ddof == 0)
			throw std::invalid_argument("Size minus degree of freedom is 0.");

		ContType<double, std::allocator<double>> x_centered = center(x);
		double sxx = std::inner_product(x_centered.begin(), x_centered.end(), x_centered.begin(), 0.0);
		double std = std::sqrt(sxx / static_cast<double>(x.size() - ddof));

		ContType<double, std::allocator<double>> z(size);
		std::transform(x_centered.begin(), x_centered.end(), z.begin(), std::bind2nd(std::divides<double>(), std));
		return z;
	}

	// --- CorrelationFunctions --- //

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
		if (x.size() != y.size())
			throw std::invalid_argument("Inputs have not the same size.");
		if (x.size() == 0)
			throw std::invalid_argument("Inputs have not enough values for pearsonr.");

		ContType<double, std::allocator<double>> x_centered = Stats::center(x);
		ContType<double, std::allocator<double>> y_centered = Stats::center(y);

		double sxx = std::inner_product(x_centered.begin(), x_centered.end(), x_centered.begin(), 0.0);
		double syy = std::inner_product(y_centered.begin(), y_centered.end(), y_centered.begin(), 0.0);
		double sxy = std::inner_product(x_centered.begin(), x_centered.end(), y_centered.begin(), 0.0);

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
		if (x.size() != y.size())
			throw std::invalid_argument("Inputs have not the same size.");
		if (x.size() == 0)
			throw std::invalid_argument("Inputs have not enough values for spearmanr.");

		ContType<unsigned int, std::allocator<unsigned int>> xr = Stats::rankdata(x);
		ContType<unsigned int, std::allocator<unsigned int>> yr = Stats::rankdata(y);

		double r = pearsonr(xr, yr);
		return r;
	}

	// --- StatisticalTests --- //

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