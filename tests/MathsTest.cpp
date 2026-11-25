#include "CppUnitTest.h"
#include <list>
#include <vector>
#include <deque>
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

		TEST_METHOD(IsPositiveTest)
		{
			std::vector<int> x = { 1, 2, 3 };
			Assert::IsTrue(Maths::is_positive(x));

			x = { 0, 1 };
			Assert::IsFalse(Maths::is_positive(x));
		}

		TEST_METHOD(ReciprocalTest)
		{
			std::vector<float> x = { -1.0f, 1.0f, 2.0f };
			std::vector<double> s = Maths::reciprocal(x);
			for (size_t i = 0; i < x.size(); ++i)
				Assert::AreEqual(1. / x[i], s[i]);
		}
		TEST_METHOD(PowerTest)
		{
			std::vector<int> x = { 0, 1, 2 };
			double p = 3.2;
			std::vector<double> s = Maths::power(x, p);
			for (size_t i = 0; i < x.size(); ++i)
				Assert::AreEqual(std::pow(x[i], p), s[i]);
		}
		TEST_METHOD(LogTest)
		{
			std::vector<float> x = { 1.0f, 2.3f, 18.4f };
			std::vector<double> s = Maths::log(x);
			for (size_t i = 0; i < x.size(); ++i)
				Assert::AreEqual(std::log(static_cast<double>(x[i])), s[i]);
		}
		TEST_METHOD(ExpTest)
		{
			std::vector<float> x = { 0.0f, 1.0f, -2.0f };
			std::vector<double> s = Maths::exp(x);
			for (size_t i = 0; i < x.size(); ++i)
				Assert::AreEqual(std::exp(static_cast<double>(x[i])), s[i]);
		}

		TEST_METHOD(SetTest)
		{
			std::vector<int> x = { 0, 1, 0, 0, 1, 1, 1, 1, 1, 1 };
			std::vector<int> s = Maths::set(x);
			Assert::AreEqual<size_t>(2, s.size());
			Assert::AreEqual(0, s[0]);
			Assert::AreEqual(1, s[1]);
		}

		TEST_METHOD(RankdataOfVectorOfDoublesTest)
		{
			std::vector<double> x = { 0.0, 1.0, -1.0, 2.0 };
			std::vector<unsigned int> r_act = Maths::rankdata(x);
			std::vector<unsigned int> r_exp = { 2, 0, 1, 3 };
			for (size_t i = 0; i < r_exp.size(); ++i)
				Assert::AreEqual(r_exp[i], r_act[i]);
		}
		TEST_METHOD(RankdataOfListOfIntegersTest)
		{
			std::list<int> x = { 0, 1, -1, 2 };
			std::list<unsigned int> r_act = Maths::rankdata(x);
			std::list<unsigned int> r_exp = { 2, 0, 1, 3 };
			for (auto rit = r_exp.begin(), rfit = r_act.begin(); rit != r_exp.end() && rfit != r_act.end(); ++rit, ++rfit)
				Assert::AreEqual(*rit, *rfit);
		}
	};
}
