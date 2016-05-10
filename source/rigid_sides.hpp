#ifndef RIGID_SIDES_HPP
#define RIGID_SIDES_HPP 1

#include "source/newtonian/two_dimensional/GhostPointGenerator.hpp"

/*! \brief Class for creating computationalcells of ghost points for rigid walls
\author Elad Steinberg
*/
class RigidSides : public GhostPointGenerator
{
public:
	boost::container::flat_map<size_t, ComputationalCell> operator() (const Tessellation& tess, const vector<ComputationalCell>& cells,
		double time, TracerStickerNames const& tracerstickernames) const;

	Slope GetGhostGradient(const Tessellation& tess,
		const vector<ComputationalCell>& cells, const vector<Slope>& gradients,
		size_t ghost_index, double time, const Edge& edge, TracerStickerNames const& tracerstickernames) const;
};

#endif // RIGID_SIDES_HPP
