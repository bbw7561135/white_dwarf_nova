#include "source/newtonian/two_dimensional/time_step_function.hpp"
#include "nuclear_burn.hpp"

class BurnTime: public TimeStepFunction
{
public:

  BurnTime(const NuclearBurn& nb,
	   const double cfl);

  double operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const vector<Vector2D>& point_velocities,
   const double time,
   TracerStickerNames const& tracerstickernames) const;

private:
  const NuclearBurn& nb_;
  const double cfl_;
};
