#include "CppUnitTest.h"
#include <vector>
#include <deque>
#include "MovingZscoreDetector.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;


namespace StatsSimpleCppTests
{
	TEST_CLASS(MovingZscoreDetectorTests)
	{
	public:

		double zscore_threshold = 2.5;
		int sample_init = 3;

		TEST_METHOD(MovingZscoreDetectorFitDetectOnVectorOfFloatsTest)
		{
			std::vector<float> x = { -1.5f, 1.5f - 1.0f, 1.0f, -0.5f, 0.5f, 0.0f };

			MovingZscoreDetector mzd = MovingZscoreDetector(zscore_threshold, sample_init);

			for (size_t i = 0; i < x.size(); i++)
				mzd.partial_fit(x[i], 0.2);

			Assert::AreEqual(1, mzd.detect(5.0));
			Assert::AreEqual(0, mzd.detect(0.0));
		}

		TEST_METHOD(MovingZscoreDetectorFitDetectOnDequeOfIntegersTest)
		{
			std::deque<int> x = { -2, 2, -1, 1, 0, 0, 1 };

			MovingZscoreDetector mzd = MovingZscoreDetector(zscore_threshold, sample_init);

			for (size_t i = 0; i < x.size(); i++)
				mzd.partial_fit(x[i], 0.5);

			Assert::AreEqual(1, mzd.detect(5.0));
			Assert::AreEqual(0, mzd.detect(0.0));
		}
	};
}
