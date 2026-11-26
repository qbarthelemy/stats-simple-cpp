#include "CppUnitTest.h"
#include <list>
#include <vector>
#include <deque>
#include "SimpleLinearRegression.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;


namespace StatsSimpleCppTests
{
	TEST_CLASS(SimpleLinearRegressionTests)
	{
	public:

		TEST_METHOD(SimpleLinearRegressionOnListOfDoublesTest)
		{
			std::list<double> x = { 55.1, 13, -75, 17.5, 89 }, y = { -2, 15, 7.2, 10, -1 };
			SimpleLinearRegression slr = SimpleLinearRegression();
			slr.fit(x, y);
			x = { 42 };
			y = slr.predict(x);
		}

		TEST_METHOD(SimpleLinearRegressionFitPredictOnVectorOfIntegersTest)
		{
			std::vector<int> x = { -1, 0, 1, 2, 3, 4 };
			std::vector<int> y = { 42, 41, 40, 39, 38, 37 };
			SimpleLinearRegression slr = SimpleLinearRegression();
			slr.fit(x, y);
			x = { 5 };
			std::vector<double> z = slr.predict(x);
			Assert::AreEqual(36.0, z[0]);
		}

		TEST_METHOD(SimpleLinearRegressionFitScoreOnDequeOfDoublesTest)
		{
			std::deque<double> x = { -1.0, 0.0, 1.0, 2.0, 3.0, 4.0 };
			std::deque<double> y = { 42.0, 41.0, 40.0, 39.0, 38.0, 37.0 };
			SimpleLinearRegression slr = SimpleLinearRegression();
			slr.fit(x, y);
			double r2 = slr.score(x, y);
			Assert::AreEqual(1.0, r2);
		}
	};
}
