#include <cmath>
#include <numeric>
#include <functional>
#include <algorithm>


namespace Maths
{
	// --- Arithmetic --- //

	/**
	* Greatest common divisor.
	*
	* @tparam ValType The numeric data type of the inputs.
	*
	* @param m, n Input values.
	*
	* @return Greatest common divisor of `m` and `n`.
	*/
	template<typename ValType>
	ValType gcd(ValType m, ValType n)
	{
		m = std::abs(m);
		n = std::abs(n);

		ValType t;
		if (m > n) // ensure n > m
		{
			t = n;
			n = m;
			m = t;
		}
		while (n != 0)
		{
			t = m % n;
			m = n;
			n = t;
		}

		return m;
	}

	/**
	* Factorial.
	*
	* @tparam ValType The numeric data type of the input.
	*
	* @param n Input value.
	*
	* @return Factorial of `n`.
	*/
	template<typename ValType>
	ValType factorial(ValType n)
	{
		n = std::abs(n);

		ValType fact = 1;
		for (ValType c = 1; c <= n; ++c)
			fact *= c;
		return fact;
	}

	// --- BasicFunctions --- //

	/**
	* Product of values.
	*
	* @tparam ContType The type of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Product of values of `x`, whose data type is the same as `x`.
	*/
	template<typename ContType>
	typename ContType::value_type prod(const ContType& x)
	{
		if (x.size() == 0)
			throw std::invalid_argument("Input has not enough values for prod.");

		ContType::value_type p = std::accumulate(x.begin(), x.end(), static_cast<ContType::value_type>(1), std::multiplies<ContType::value_type>());
		return p;
	}

	/**
	* Absolute value, element-wise.
	*
	* @tparam ContType The type of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Output sequence container containing the element-wise absolute value of `x`.
	*/
	template<typename ContType>
	ContType absolute(const ContType& x)
	{
		ContType x_abs(x.size());
		std::transform(x.begin(), x.end(), x_abs.begin(), static_cast<ContType::value_type(*)(ContType::value_type)>(&std::abs));
		return x_abs;
	}

	/**
	* Inverse, element-wise.
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Output sequence container containing the element-wise inverse of `x`,
	* whose data type is double to keep the maximum of numerical precision.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	typename ContType<double, std::allocator<double>> inverse(const ContType<ValType, Alloc>& x)
	{
		//TODO: check non zero

		ContType<double, std::allocator<double>> x_inverse(x.size());
		std::transform(x.begin(), x.end(), x_inverse.begin(), [](const ValType e) { return 1.0 / static_cast<double>(e); });
		return x_inverse;
	}

}