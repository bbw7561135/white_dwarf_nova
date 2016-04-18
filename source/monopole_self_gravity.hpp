#ifndef MONOPOLE_SELF_GRAVITY_HPP
#define MONOPOLE_SELF_GRAVITY_HPP 1

#include "source/newtonian/two_dimensional/SourceTerm.hpp"

class MonopoleSelfGravity: public SourceTerm
{
public:

  MonopoleSelfGravity(const vector<double>& sample_radii,
		      double gravitation_constant,
		      const pair<double,double>& section_angles);

  vector<Extensive> operator()
  (const Tessellation& tess,
   const PhysicalGeometry& /*pg*/,
   const CacheData& cd,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& /*fluxes*/,
   const vector<Vector2D>& /*point_velocities*/,
   const double /*t*/,
   const TracerStickerNames& tsn) const;

private:
  const vector<double> sample_radii_;
  const double gravitation_constant_;
  const double section2shell_;
};

#endif // MONOPOLE_SELF_GRAVITY_HPP
