#include "CppUnitTest.h"
#include <cmath>
#include <list>
#include <vector>
#include <deque>
#include "Stats.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;


namespace StatsSimpleCppTests
{
	TEST_CLASS(StatsTests)
	{
	public:

		TEST_METHOD(MeanOfVectorOfShortsTest)
		{
			std::vector<short> x = { 1, 2, 3 };
			Assert::AreEqual(2.0, stats::mean(x));
		}
		TEST_METHOD(MeanOfDequeOfFloatsTest)
		{
			std::deque<float> x = { -1.0f, 12.0f, 1.0f };
			Assert::AreEqual(4.0, stats::mean(x), 1e-16);
		}
		TEST_METHOD(MeanOfListOfDoublesTest)
		{
			std::list<double> x = { -1.0, 3e-10, 1.0 };
			Assert::AreEqual(1e-10, stats::mean(x), 1e-16);
		}
		TEST_METHOD(MeanOfVoidVectorOfDoublesTest)
		{
			Assert::ExpectException<std::invalid_argument>([]
				{
					std::vector<double> c;
					stats::mean(c);
				});
		}

		TEST_METHOD(CenteringOfVectorOfIntegersTest)
		{
			std::vector<int> x = { 0, 1, 2, 3, 4 };
			std::vector<double> center = stats::center(x);
			Assert::AreEqual(0.0, stats::mean(center), 1e-10);
		}
		TEST_METHOD(CenteringOfDequeOfFloatsTest)
		{
			std::deque<float> x = { -3.57f, 3.42f, 17.99f };
			std::deque<double> center = stats::center(x);
			Assert::AreEqual(0.0, stats::mean(center), 1e-10);
		}
		TEST_METHOD(CenteringOfListOfDoublesTest)
		{
			std::list<double> x = { -1.1, 1.2, 2.17, 3.42, 17.99 };
			std::list<double> center = stats::center(x);
			Assert::AreEqual(0.0, stats::mean(center), 1e-10);
		}

		TEST_METHOD(HMeanOfDequeOfFloatsTest)
		{
			std::deque<float> x = { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f };
			Assert::AreEqual(2.6997245179063363, stats::hmean(x), 1e-3);
		}
		TEST_METHOD(GMeanOfListOfDoublesTest)
		{
			std::list<double> x = { 1, 2, 3, 4, 5, 6, 7 };
			Assert::AreEqual(3.3800151591412964, stats::gmean(x), 1e-6);
		}
		TEST_METHOD(PMeanOfVectorOfLongLongsTest)
		{
			std::vector<long long> x = { 1, 2, 3, 4, 5, 6, 7 };
			Assert::AreEqual(4.157111214492084, stats::pmean(x, 1.3), 1e-6);
		}

		TEST_METHOD(VarianceOfVectorOfFloatsTest)
		{
			std::vector<float> x = { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f };
			Assert::AreEqual(2.0, stats::var(x));
		}
		TEST_METHOD(VarianceOfDequeOfDoublesTest)
		{
			std::deque<double> x = { 1.0, 2.0, 3.0, 4.0 };
			Assert::AreEqual(1.25, stats::var(x));
		}
		TEST_METHOD(VarianceOfVoidVectorOfDoublesTest)
		{
			Assert::ExpectException<std::invalid_argument>([]
				{
					std::vector<double> x;
					stats::var(x);
				});
		}

		TEST_METHOD(StdOfVectorOfDoublesTest)
		{
			std::vector<double> x = { 0.0, 1.0, 2.0, 3.0, 4.0 };
			Assert::AreEqual(std::sqrt(2.0), stats::std(x));
		}

		TEST_METHOD(GStdOfListOfDoublesTest)
		{
			std::list<double> x = { 1.82424757f, 1.22436866f, 1.13183117f };
			Assert::AreEqual(1.23241, stats::gstd(x), 1e-5);
		}

		TEST_METHOD(SkewnessOfListOfDoublesTest)
		{
			std::list<double> x = { 2, 8, 0, 4, 1, 9, 9, 0 };
			Assert::AreEqual(0.2650554122698573, stats::skew(x), 1e-5);
		}

		TEST_METHOD(KurtosisOfListOfFloatsTest)
		{
			std::list<float> x = { -4.1f, 1.2f, 2.0f, 4.0f, 13.42f, 117.99f };
			Assert::AreEqual(4.09302, stats::kurtosis(x), 1e-5);
		}

		TEST_METHOD(MedianOfVectorOfShortsTest)
		{
			std::vector<short> x = { -1, 1, 2, 3, 4, 17 };
			Assert::AreEqual(2.5, stats::median(x));
		}
		TEST_METHOD(MedianOfListOfFloatTest)
		{
			std::list<float> x = { -4.1f, 1.2f, 2.0f, 4.0f, 13.42f, 117.99f };
			Assert::AreEqual(3.0, stats::median(x));
		}
		TEST_METHOD(MedianOfDequeOfDoublesTest)
		{
			std::deque<double> x = { -1.1, 1.2, 2.17, 3.42, 17.99 };
			Assert::AreEqual(2.17, stats::median(x));
		}
		TEST_METHOD(MedianOfVoidVectorOfDoublesTest)
		{
			Assert::ExpectException<std::invalid_argument>([]
				{
					std::vector<double> x;
					stats::median(x);
				});
		}

		TEST_METHOD(MedianAbsDevOfDequeOfDoublesTest)
		{
			std::deque<double> x = { -1.1, 1.2, 2.17, 3.42, 17.99 };
			Assert::AreEqual(1.25, stats::median_abs_deviation(x));
		}

		TEST_METHOD(ZscoreOfVectorOfDoublesTest)
		{
			std::vector<double> x = { -0.7972, 0.0767, 0.4383, 0.7866, 0.8091,
				0.1954, 0.6307, 0.6599, 0.1065, 0.0508 };
			std::vector<double> z_act = stats::zscore(x);
			std::vector<double> z_exp = { -2.37471, -0.47582 , 0.309898, 1.06672, 1.11561,
				-0.217898,  0.727962, 0.791411, -0.411068, -0.532098 };
			for (size_t i = 0; i < z_exp.size(); ++i)
				Assert::AreEqual(z_exp[i], z_act[i], 1e-5);
		}

		TEST_METHOD(GZscoreOfVectorOfDoublesTest)
		{
			std::vector<double> x = { 0.7972,  0.0767,  0.4383,  0.7866,  0.8091,
				0.1954,  0.6307,  0.6599,  0.1065,  0.0508 };
			std::vector<double> z_act = stats::gzscore(x);
			std::vector<double> z_exp = { 0.932574, -1.35761, 0.347408, 0.91948, 0.947068,
				-0.442842, 0.703404, 0.747676, -1.03652, -1.76064 };
			for (size_t i = 0; i < z_exp.size(); ++i)
				Assert::AreEqual(z_exp[i], z_act[i], 1e-5);
		}

		TEST_METHOD(PearsonCorrelationPositiveTest)
		{
			std::vector<double> v1 = { -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
			std::vector<double> v2 = { 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0 };
			Assert::AreEqual(1.0, stats::pearsonr(v1, v2));
		}
		TEST_METHOD(PearsonCorrelationNegativeTest)
		{
			std::list<int> v1 = { 0, 1, 2, 3, 4, 5 };
			std::list<int> v2 = { 22, 21, 20, 19, 18, 17 };
			Assert::AreEqual(-1.0, stats::pearsonr(v1, v2));
		}
		TEST_METHOD(PearsonCorrelationDifferentSizeTest)
		{
			Assert::ExpectException<std::invalid_argument>([]
				{
					std::vector<double> v1 = { 0.0, 1.0 };
					std::vector<double> v2 = { 17.0 };
					stats::pearsonr(v1, v2);
				});
		}

		TEST_METHOD(SpearmanCorrelationPositiveTest)
		{
			std::vector<double> v1 = { 1, 2, 3, 4, 5 };
			std::vector<double> v2 = { 5, 6, 7, 8, 7 };
			Assert::AreEqual(0.8207826816681233, stats::spearmanrho(v1, v2));
		}
		TEST_METHOD(SpearmanCorrelationNegativeTest)
		{
			std::list<int> v1 = { 0, 1, 2, 3, 4, 5 };
			std::list<int> v2 = { 22, 21, 20, 19, 18, 17 };
			Assert::AreEqual(-1.0, stats::spearmanrho(v1, v2));
		}

		TEST_METHOD(AccuracyVectorOfIntegersTest)
		{
			std::vector<int> y1 = { 0, 0, 0, 1, 1, 1 };
			std::vector<int> y2 = { 0, 0, 0, 1, 1, 1 };
			Assert::AreEqual(1.0, stats::accuracy_score(y1, y2));
		}
		TEST_METHOD(AccuracyListOfIntegersTest)
		{
			std::list<int> y1 = { 0, 0, 1, 1 };
			std::list<int> y2 = { 0, 1, 1, 0 };
			Assert::AreEqual(0.5, stats::accuracy_score(y1, y2));
		}
	};
}
