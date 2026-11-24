#include "CppUnitTest.h"
#include <vector>
#include "MovingZscoreDetector.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;


namespace StatsSimpleCppTests
{
	TEST_CLASS(MovingZscoreDetectorTests)
	{
	public:

		TEST_METHOD(MovingZscoreDetectorFitDetect)
		{
			std::vector<float> x = { -1.5f, 1.5f - 1.0f, 1.0f, -0.5f, 0.5f, 0.0f };
			//std::vector<int> x = { -2, 2, -1, 1, 0, 0, 1 };

			double zscore_threshold = 2.5;
			int sample_init = static_cast<int>(x.size() / 2);

			MovingZscoreDetector mzd = MovingZscoreDetector(zscore_threshold, sample_init);

			for (size_t i = 0; i < x.size(); i++)
				mzd.partial_fit(x[i], 0.2);

			Assert::AreEqual(1, mzd.detect(5.0));
			Assert::AreEqual(0, mzd.detect(0.0));
		}
	};
}
