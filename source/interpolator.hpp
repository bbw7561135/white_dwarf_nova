#ifndef INTERPOLATOR_HPP
#define INTERPOLATOR_HPP 1

#include <vector>

using std::vector;

class Interpolator
{
public:

	Interpolator(void);

  Interpolator(const vector<double>& x_list,
	       const vector<double>& y_list);

  Interpolator& operator=(Interpolator const& other);

  double operator()(double x) const;

private:
  vector<double> x_list_;
  vector<double> y_list_;
};

#endif // INTERPOLATOR_HPP
