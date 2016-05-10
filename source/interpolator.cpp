#include "interpolator.hpp"
#include "vector_utils.hpp"
#include <cassert>

using std::size_t;

Interpolator::Interpolator(const vector<double>& x_list,
			   const vector<double>& y_list):
  x_list_(x_list), y_list_(y_list)
{
  assert(is_strictly_increasing(x_list_));
  assert(x_list_.size()==y_list_.size());
}

double Interpolator::operator()(double x) const
{
  const double x_list_front = x_list_.front();
  const double x_list_back = x_list_.back();
  assert(x>x_list_front);
  assert(x<x_list_back);       
  //      assert(x>x_list_.front());
  //      assert(x<x_list_.back());
  for(size_t i=1;i<x_list_.size();++i){
    if(x_list_.at(i)>x)
      return y_list_.at(i-1) + (y_list_.at(i)-y_list_.at(i-1))*
	(x-x_list_.at(i-1))/(x_list_.at(i)-x_list_.at(i-1));
  }
  throw "point outside bound";
}

Interpolator& Interpolator::operator=(Interpolator const& other)
{
	if (this == &other)
		return *this;
	x_list_ = other.x_list_;
	y_list_ = other.y_list_;
	return *this;
}

Interpolator::Interpolator(void) :x_list_(vector<double>()), y_list_(vector<double>()) {}