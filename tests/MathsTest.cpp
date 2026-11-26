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
			Assert::AreEqual(4, maths::gcd(28, 32));
		}

		TEST_METHOD(FactorialTest)
		{
			Assert::AreEqual(6, maths::factorial(3));
		}

		TEST_METHOD(IsPositiveTest)
		{
			std::vector<int> x = { 1, 2, 3 };
			Assert::IsTrue(maths::is_positive(x));

			x = { 0, 1 };
			Assert::IsFalse(maths::is_positive(x));
		}

		TEST_METHOD(ProdTest)
		{
			std::vector<int> x = { 1, 2, 3 };
			Assert::AreEqual(6, maths::prod(x));
		}

		TEST_METHOD(LinearTest)
		{
			std::vector<int> x = { 1, 2, 3 };
			double a = 2, b = -2;
			std::vector<double> y_act = maths::linear(x, a, b);
			std::vector<double> y_exp = { 0, 2, 4 };
			for (size_t i = 0; i < y_exp.size(); ++i)
				Assert::AreEqual(y_exp[i], y_act[i]);
		}

		TEST_METHOD(AbsoluteTest)
		{
			std::deque<float> x = { -1.0f, 2.0, -3.0f };
			std::deque<float> x_act = maths::absolute(x);
			for (size_t i = 0; i < x_act.size(); ++i)
				Assert::AreEqual(std::abs(x[i]), x_act[i]);
		}

		TEST_METHOD(ReciprocalTest)
		{
			std::vector<float> x = { -1.0f, 1.0f, 2.0f };
			std::vector<double> y = maths::reciprocal(x);
			for (size_t i = 0; i < x.size(); ++i)
				Assert::AreEqual(1. / x[i], y[i]);
		}
		TEST_METHOD(PowerTest)
		{
			std::vector<int> x = { 0, 1, 2 };
			double p = 3.2;
			std::vector<double> y = maths::power(x, p);
			for (size_t i = 0; i < x.size(); ++i)
				Assert::AreEqual(std::pow(x[i], p), y[i]);
		}
		TEST_METHOD(LogTest)
		{
			std::vector<float> x = { 1.0f, 2.3f, 18.4f };
			std::vector<double> y = maths::log(x);
			for (size_t i = 0; i < x.size(); ++i)
				Assert::AreEqual(std::log(static_cast<double>(x[i])), y[i]);
		}
		TEST_METHOD(ExpTest)
		{
			std::vector<float> x = { 0.0f, 1.0f, -2.0f };
			std::vector<double> y = maths::exp(x);
			for (size_t i = 0; i < x.size(); ++i)
				Assert::AreEqual(std::exp(static_cast<double>(x[i])), y[i]);
		}

		TEST_METHOD(SigmoidTest)
		{
			std::deque<float> x = { -1.0f, 0.0f, 1.0f };
			std::deque<double> y_act = maths::sigmoid(x);
			std::deque<double> y_exp = { 0.268941, 0.5, 0.731059 };
			for (size_t i = 0; i < x.size(); ++i)
				Assert::AreEqual(y_exp[i], y_act[i], 1e-5);
		}

		TEST_METHOD(SetTest)
		{
			std::vector<int> x = { 0, 1, 0, 0, 1, 1, 1, 1, 1, 1 };
			std::vector<int> y = maths::set(x);
			Assert::AreEqual<size_t>(2, y.size());
			Assert::AreEqual(0, y[0]);
			Assert::AreEqual(1, y[1]);
		}

		TEST_METHOD(RankdataOfVectorOfDoublesTest)
		{
			std::vector<double> x = { 0.0, 1.0, -1.0, 2.0 };
			std::vector<unsigned int> r_act = maths::rankdata(x);
			std::vector<unsigned int> r_exp = { 2, 0, 1, 3 };
			for (size_t i = 0; i < r_exp.size(); ++i)
				Assert::AreEqual(r_exp[i], r_act[i]);
		}
		TEST_METHOD(RankdataOfListOfIntegersTest)
		{
			std::list<int> x = { 0, 1, -1, 2 };
			std::list<unsigned int> r_act = maths::rankdata(x);
			std::list<unsigned int> r_exp = { 2, 0, 1, 3 };
			for (auto rit = r_exp.begin(), rfit = r_act.begin(); rit != r_exp.end() && rfit != r_act.end(); ++rit, ++rfit)
				Assert::AreEqual(*rit, *rfit);
		}
	};
}
