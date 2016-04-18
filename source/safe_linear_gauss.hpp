/*! \file LinearGaussImproved.hpp
\brief Linear interpolation that guarantees compliance with the equation of state and calcualtes the GG gradient from the CM of the cell
\author Elad Steinberg
*/

#ifndef SAFE_LINEAR_GAUSS_HPP
#define SAFE_LINEAR_GAUSS_HPP 1

#include "fermi_table.hpp"
#include "source/newtonian/two_dimensional/spatial_reconstruction.hpp"
#include <cmath>
#include "source/newtonian/two_dimensional/GhostPointGenerator.hpp"

//! \brief Linear gauss interpolation
class SafeLinearGauss : public SpatialReconstruction
{
public:

	/*! \brief Class constructor
	\param eos Equation of state
	\param slf Slope limiter flag
	\param delta_v The GradV*L/Cs ratio needed for slope limiter
	\param theta The theta from tess in slope limiter.
	\param delta_P The pressure ratio for shock detection
	\param ghost The ghost point generator
	\param flat_tracers Names of tracers for which the slope is always zero
	*/
	SafeLinearGauss
		(FermiTable const& eos,
		GhostPointGenerator const& ghost,
		bool slf = true,
		double delta_v = 0.2,
		double theta = 0.5,
		double delta_P = 0.7,
		const vector<string>& flat_tracers =
		vector<string>());

	void operator() 
	(const Tessellation& tess,
	 const vector<ComputationalCell>& cells,
	 double time,
	 vector<pair<ComputationalCell, ComputationalCell> >& res,
	 const TracerStickerNames& tsn) const;

	/*! \brief Interpolates a cell
	\param cell The primitives of the cell
	\param cell_index The index of the cell
	\param cm The cell's center of mass
	\param target The location of the interpolation
	\return The interpolated value
	*/
	ComputationalCell Interp
	(const ComputationalCell& cell, 
	 const Slope& slope,
	 const Vector2D& target,
	 const Vector2D& cm) const;

	/*!
	\brief Returns the gradients
	\return The gradients
	*/
  vector<Slope>& GetSlopes(void) const;

	/*!
	\brief Returns the unsloped limtied gradients
	\return The gradients
	*/
	vector<Slope>& GetSlopesUnlimited(void)const;
private:
  const FermiTable& eos_;
	GhostPointGenerator const& ghost_;
	mutable vector<Slope> rslopes_;
	mutable vector<Slope> naive_rslopes_;
	const bool slf_;
	const double shockratio_;
	const double diffusecoeff_;
	const double pressure_ratio_;
	const vector<string> flat_tracers_;

	SafeLinearGauss
		(const SafeLinearGauss& origin);

	SafeLinearGauss& operator=
		(const SafeLinearGauss& origin);

};
#endif // SAFE_LINEAR_GAUSS_HPP
