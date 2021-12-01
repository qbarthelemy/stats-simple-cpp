#include <cmath>


namespace Maths
{

	/**
	* Greatest common divisor.
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
	*/
	template<typename ValType>
	ValType factorial(ValType n)
	{
		ValType f = 1;

		for (ValType c = 1; c <= n; ++c)
			f = f * c;

		return f;
	}

}