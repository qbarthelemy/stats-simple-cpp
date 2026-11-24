#include "CppUnitTest.h"
#include <vector>
#include "Maths.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;


namespace StatsSimpleCppTests
{
	TEST_CLASS(MathsTests)
	{
	public:

		TEST_METHOD(GcdTest)
		{
			Assert::AreEqual(4, Maths::gcd(28, 32));
		}

		TEST_METHOD(FactorialTest)
		{
			Assert::AreEqual(6, Maths::factorial(3));
		}

		TEST_METHOD(PowerTest)
		{
			std::vector<int> x = { 0, 1, 0, 0, 1, 1, 1, 1, 1, 1 };
			std::vector<double> s = Maths::power(x, 3.2);
		}
		TEST_METHOD(LogTest)
		{
			std::vector<float> x = { 1.0f, 2.3f, 18.4f };
			std::vector<double> s = Maths::log(x);
		}
		TEST_METHOD(ExpTest)
		{
			std::vector<float> x = { 1.0f, 2.3f, 18.4f };
			std::vector<double> s = Maths::exp(x);
		}

		TEST_METHOD(SetTest)
		{
			std::vector<int> x = { 0, 1, 0, 0, 1, 1, 1, 1, 1, 1 };
			std::vector<int> s = Maths::set(x);
			Assert::AreEqual<size_t>(2, s.size());
			Assert::AreEqual(0, s[0]);
			Assert::AreEqual(1, s[1]);
		}
	};
}
