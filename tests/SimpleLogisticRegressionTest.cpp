#include "CppUnitTest.h"
#include <list>
#include <vector>
#include <deque>
#include "SimpleLogisticRegression.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;


namespace StatsSimpleCppTests
{
	TEST_CLASS(SimpleLogisticRegressionTests)
	{
	public:

		TEST_METHOD(SimpleLogisticRegressionFitPredictOnVectorOfIntegersTest)
		{
			std::vector<int> x = { -1, 0, 1, 2, 3, 4 };
			std::vector<int> y = { 0, 0, 0, 1, 1, 1 };
			SimpleLogisticRegression slr = SimpleLogisticRegression(0.1, 0.1, 5000);
			slr.fit(x, y);
			x = { 5 };
			std::vector<double> proba = slr.predict_proba(x);
			Assert::AreEqual(0.992221, proba[0], 1e-5);
			std::vector<int> pred = slr.predict(x);
			Assert::AreEqual(1, pred[0]);
		}
		TEST_METHOD(SimpleLogisticRegressionFitScoreOnDequeOfFloatsTest)
		{
			std::deque<float> x = { -1.0f, 0.0f, 1.0f, 12.0f, 13.0f, 14.0f };
			std::deque<int> y = { 1, 1, 1, 0, 0, 0 };
			SimpleLogisticRegression slr = SimpleLogisticRegression(0.1, 0.01, 5000);
			slr.fit(x, y);
			double acc = slr.score(x, y);
			Assert::AreEqual(1.0, acc);
		}
		TEST_METHOD(SimpleLogisticRegressionFitScoreOnListOfDoublesTest)
		{
			std::list<double> x = { -10.2, 0.1, 1.0, 12.0, 300.0, 140.0 };
			std::list<int> y = { 1, 1, 1, 0, 0, 0 };
			SimpleLogisticRegression slr = SimpleLogisticRegression(0.1, 0.01, 5000);
			slr.fit(x, y);
			double acc = slr.score(x, y);
			Assert::AreEqual(1.0, acc);
		}
	};
}
