#ifndef CORE_ATMOSPHERE_GRAVITY_HPP
#define CORE_ATMOSPHERE_GRAVITY_HPP 1

#include "source/newtonian/two_dimensional/source_terms/ConservativeForce.hpp"
#include "interpolator.hpp"

using std::vector;

class EnclosedMassCalculator
{
public:
	double operator()(double radius) const;
	
	void ReCalcEnclosedMass(const double core_mass,const vector<pair<double,double> >& mass_radius_list,
		const vector<double>& sample_radii,const double section2shell);

private:
	Interpolator interpolator_;
};


class AccelerationCalculator : public Acceleration
{
public:

	AccelerationCalculator(const double gravitation_constant,const EnclosedMassCalculator& emc);

	Vector2D operator()
	(const Tessellation& tess,
	const vector<ComputationalCell>& cells,
	const vector<Extensive>& fluxes,
	const double time,
	const int point,
	TracerStickerNames const& tracerstickernames) const;

private:
	const double gravitation_constant_;
	const EnclosedMassCalculator& emc_;
};

class CoreAtmosphereGravity: public SourceTerm
{
public:

	CoreAtmosphereGravity
	(const double core_mass,
	const vector<double>& sample_radii,
	const double gravitation_constant,
	const pair<double,double>& sector_angles);

	vector<Extensive> operator()
	(const Tessellation& tess,
	const PhysicalGeometry& pg,
	const CacheData& cd,
	const vector<ComputationalCell>& cells,
	const vector<Extensive>& fluxes,
	const vector<Vector2D>& point_velocities,
	const double time,
	const TracerStickerNames& tsn) const;

	double getCoreMass(void) const;

	const vector<double>& getSampleRadii(void) const;

	double getSection2Shell(void) const;

	double getGravitationConstant(void) const;

private:
	const double core_mass_;
	const vector<double> sample_radii_;
	const double gravitation_constant_;
	const double section2shell_;
	mutable EnclosedMassCalculator emc_;
	AccelerationCalculator acc_;
};

#endif // CORE_ATMOSPHERE_GRAVITY_HPP
