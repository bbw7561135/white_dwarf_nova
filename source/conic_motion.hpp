#ifndef CONIC_MOTION_HPP
#define CONIC_MOTION_HPP 1

#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "circular_section.hpp"
class ConicMotion: public PointMotion
{
public:

	ConicMotion(const CircularSection& domain,const EquationOfState& eos);

	vector<Vector2D> operator()(const Tessellation& tess,const vector<ComputationalCell>& cells,double time,
	TracerStickerNames const& tracerstickernames) const;

	vector<Vector2D> ApplyFix(Tessellation const& tess, vector<ComputationalCell> const& cells, double time,
	double dt, vector<Vector2D> const& velocities, TracerStickerNames const& tracerstickernames)const;
private:
	CircularSection const& circ_section_;
	const Lagrangian bpm_;
	const RoundCells pm_;
};

#endif // CONIC_MOTION_HPP
