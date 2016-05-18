#ifndef BURN_TIME_HPP
#define BURN_TIME_HPP 1

#include "source/newtonian/two_dimensional/time_step_function.hpp"
#include "nuclear_burn.hpp"

class BurnTime: public TimeStepFunction
{
public:

  BurnTime(const double cfl);

  void attach(NuclearBurn* nb) const
  {
    nb_ = nb;
  }

  double operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const vector<Vector2D>& point_velocities,
   const double time,
   TracerStickerNames const& tracerstickernames) const;

private:
  mutable NuclearBurn* nb_;
  const double cfl_;
};

#endif // BURN_TIME_HPP
