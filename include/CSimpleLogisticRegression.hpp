#include <stdexcept>
#include <numeric>
#include <algorithm>
#include "Stats.hpp"


/**
* Simple logistic regression, for binary classification, using gradient descent.
*
* @see [sklearn.linear_model.LogisticRegression](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html)
*/
class SimpleLogisticRegression
{
public:

	/**
	* Create logistic model.
	*
	* @param learning_rate Learning rate for gradient descent.
	* @param gradient_threshold Threshold in percentage of convergence of gradient descent.
	* @param iteration_threshold Threshold on maximal number of iterations for gradient descent.
	*/
	SimpleLogisticRegression(double learning_rate = 0.001, double gradient_threshold = 0.01, int iteration_threshold = 100)
	{ 
		_coeff = 1E42, _intercept = 1E42;
		_learning_rate = learning_rate;
		_gradient_threshold = gradient_threshold;
		_iteration_threshold = iteration_threshold;
	}

	double get_coeff() { return _coeff; }

	double get_intercept() { return _intercept; }

	/**
	* Fit logistic model, using gradient descent on binary-cross entropy.
	*
	* @tparam ContType The type of the sequence containers.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container containing training values.
	* @param y Input sequence container containing binary target values.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	void fit(const ContType<ValType, Alloc>& x, const ContType<int, std::allocator<int>>& y)
	{
		if (_learning_rate <= 0)
			throw std::invalid_argument("Parameter learning_rate must be positive.");
		if (_gradient_threshold <= 0 || _gradient_threshold >=1)
			throw std::invalid_argument("Parameter gradient_threshold must be a percentage in (0, 1).");
		if (_iteration_threshold <= 0)
			throw std::invalid_argument("Parameter iteration_threshold must be positive.");

		size_t size = x.size();
		if (size != y.size())
			throw std::invalid_argument("Inputs have not the same size.");
		if (size < 2)
			throw std::invalid_argument("Inputs have not enough values for fit.");
		ContType<int, std::allocator<int>> y_set = Maths::set(y);
		if (y_set.size() != 2)
			throw std::invalid_argument("Targets must contain two classes of values.");
		if ((y_set.front() != 0 && y_set.front() != 1) || (y_set.back() != 0 && y_set.back() != 1))
			throw std::invalid_argument("Targets must contain binary values, 0 or 1.");

		double size_inv = 1.0 / size;
		ContType<int, std::allocator<int>> y_predict(size);
		ContType<double, std::allocator<double>> res(size);
		double d_coeff, d_intercept;
		int iteration_count = 0;

		_coeff = 0, _intercept = 0;
		do
		{
			y_predict = predict(x);
			std::transform(y_predict.begin(), y_predict.end(), y.begin(), res.begin(), std::minus<double>());

			d_coeff = size_inv * std::inner_product(x.begin(), x.end(), res.begin(), 0.0);
			d_intercept = size_inv * std::accumulate(res.begin(), res.end(), 0.0);
			_coeff -= _learning_rate * d_coeff;
			_intercept -= _learning_rate * d_intercept;

			iteration_count++;
		}
		while ((std::fabs(d_coeff / _coeff) > _gradient_threshold
			&& std::fabs(d_intercept / _intercept) > _gradient_threshold)
			|| iteration_count < _iteration_threshold);

		return;
	}

	/**
	* Predict using the logistic model.
	*
	* @tparam ContType The type of the sequence container.
	* @tparam ValType The numeric data type of the values of the sequence container.
	*
	* @param x Input sequence container.
	*
	* @return Output sequence container containing predicted values from `x`,
	* whose data type is int to provide a binary target value.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	ContType<int, std::allocator<int>> predict(const ContType<ValType, Alloc>& x)
	{
		ContType<double, std::allocator<double>> y_lin = Maths::linear(x, get_coeff(), get_intercept());
		ContType<double, std::allocator<double>> y_sig = Maths::sigmoid(y_lin);
		
		ContType<int, std::allocator<int>> y(x.size());
		std::transform(y_sig.begin(), y_sig.end(), y.begin(), [](double e) { return e >= 0.5 ? 1 : 0; });

		return y;
	}

	/**
	* Return the accuracy of the prediction.
	*
	* @tparam ContType The type of the sequence containers.
	* @tparam ValType The numeric data type of the values of the sequence containers.
	*
	* @param x Input sequence container containing test values.
	* @param y Input sequence container containing true binary predictions.
	*
	* @return Accuracy of `predict(x)` wrt. `y`.
	*/
	template<template<typename, typename> class ContType, typename ValType, typename Alloc>
	double score(const ContType<ValType, Alloc>& x, const ContType<int, std::allocator<int>>& y)
	{
		ContType<int, std::allocator<int>> y_predict = predict(x);
		return Stats::accuracy_score(y, y_predict);
	}

protected:

	double _coeff, _intercept;
	double _learning_rate;
	double _gradient_threshold;
	int _iteration_threshold;
};
