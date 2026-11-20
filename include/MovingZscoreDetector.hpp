#include <stdexcept>
#include <algorithm>


/**
* Moving z-score detector, to detect anomaly in an online setup.
*
* Anomaly is detected by z-score thresholding, ie when (x - mean) / std = zscore(x) > zscore_threshold
*
* @see [The Riemannian Potato Field: A Tool for Online Signal Quality Index of EEG](https://hal.science/hal-02015909/document)
* @see [Unsupervised and Adaptive Perimeter Intrusion Detector](https://hal.science/hal-03800865/document)
*/
class MovingZscoreDetector
{
public:

	/**
	* Create moving z-score detector.
	*
	* @param zscore_threshold Z-score threshold to detect anomaly.
	* @param sample_init Number of values used for initialization. No detection during initialization. 
	* @param var_min Minimal value for variance, strictly positive.
	* @param diff_min Minimal value for difference x - mean, positive. Avoid spurious detections on small raw differences.
	*/
	MovingZscoreDetector(double zscore_threshold = 2.5, int sample_init = 10, double var_min = 1E-6, double diff_min = 0)
	{
		_zscore_threshold = zscore_threshold;
		_sample_init = sample_init;
		_var_min = (var_min > 0) ? var_min : 1E-6;
		_diff_min = (diff_min >= 0) ? diff_min : 0;

		_mean = 0;
		_var = _var_min;
		_sample_count = 0;
	}

	double get_count() const { return _sample_count; }

	double get_mean() const { return _mean; }

	double get_var() const { return _var; }

	template<typename ValType>
	double get_zscore(ValType x) const { return (x - _mean) / std::sqrt(_var); }

	/**
	* Reset moving z-score detector.
	*/
	void reset()
	{
		_mean = 0;
		_var = _var_min;
		_sample_count = 0;
	}

	/**
	* Fit moving z-score detector with online update.
	*
	* @tparam ValType The numeric data type of the input.
	*
	* @param x Input value.
	* @param update_factor Update factor in (0, 1) such that:
	*   mean(t) = (1-a) * mean(t-1) + a * x(t)
	*   var(t)  = (1-a) * var(t-1)  + a * (x(t) - mean(t))^2
	* If update factor is constant, it gives an exponential moving average;
	* if it is equal to 1/t, it gives a cumulative moving average.
	*/
	template<typename ValType>
	void partial_fit(ValType x, double update_factor = 0.1)
	{
		// online initialization: cumulate value and squared value
		if (_sample_count < _sample_init)
		{
			_sample_count++;
			_mean += x;
			_var += x * x;

			// end of initialization: compute exact mean and variance (using Koenig theorem)
			if (_sample_count == _sample_init)
			{
				_mean = _mean / _sample_init;
				_var = _var / _sample_init - _mean * _mean;
			}
		}
		// online update by moving average
		else
		{
			if (update_factor <= 0 || update_factor >= 1)
				throw std::invalid_argument("Parameter update_factor must be a in (0, 1).");

			_mean = (1.0 - update_factor) * _mean + update_factor * x;
			_var = (1.0 - update_factor) * _var + update_factor * (x - _mean) * (x - _mean);

			// allow a safe use of std::sqrt() in the following
			_var = max(_var, _var_min);
		}
	}

	/**
	* Detect anomaly using moving z-score detector.
	*
	* @tparam ValType The numeric data type of the input.
	*
	* @param x Input value.
	*
	* @return Prediction of anomaly: 1 if anomaly detected, 0 otherwise.
	*/
	template<typename ValType>
	int detect(ValType x) const
	{
		if (_sample_count < _sample_init)
			return 0;

		// zscore thresholding: (x - mean) / std = z > z_th => x - mean > z_th * sqrt(var)
		double threshold = _zscore_threshold * std::sqrt(_var);
		// set a minimal value for threshold applied on x - mean
		threshold = max(threshold, _diff_min);

		if (x - _mean > threshold)
			return 1;
		else
			return 0;
	}


protected:

	double _zscore_threshold, _var_min, _diff_min, _mean, _var;
	int _sample_init, _sample_count;

};
