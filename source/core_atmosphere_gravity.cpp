#include "core_atmosphere_gravity.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#endif //RICH_MPI

namespace 
{
	vector<pair<double,double> > calc_mass_radius_list
	(const Tessellation& tess,
	const vector<ComputationalCell>& cells,
	const TracerStickerNames& tsn,
	const CacheData& cd)
	{
		vector<pair<double, double> > res;
		size_t Nloop=static_cast<size_t>(tess.GetPointNo());
		for(size_t i=0;i<Nloop;++i){
			if(safe_retrieve
					(cells.at(i).stickers,
						tsn.sticker_names,
						string("ghost")))
			continue;
			const double radius = abs(tess.GetCellCM(static_cast<int>(i)));
			const double mass = cd.volumes[i]*cells[i].density;
			res.push_back(pair<double,double>(radius,mass));
		}
		return res;
	}

	vector<double> calc_mass_in_shells
	(const vector<pair<double,double> >& mass_radius_list,
	const vector<double>& sample_points)
	{
		vector<double> res(sample_points.size(),0);
		for(vector<pair<double,double> >::const_iterator it=
		mass_radius_list.begin();
		it!=mass_radius_list.end();
		++it){
			for(size_t i=0;i<sample_points.size();++i){
				if(sample_points[i]>it->first){
					res[i] += it->second;
				}
			}
		}
#ifdef RICH_MPI
		vector<double> mpires(res.size());
		MPI_Allreduce(&res[0], &mpires[0],static_cast<int>(mpires.size()),MPI_DOUBLE,
		MPI_SUM,MPI_COMM_WORLD);
		res=mpires;
#endif //RICH_MPI
		return res;
	}

	vector<double> mult_all(double s,
	const vector<double>& v)
	{
		vector<double> res = v;
		for(size_t i=0;i<v.size();++i)
		res[i] = v[i]*s;
		return res;
	}

	vector<double> add_all(double s,
	const vector<double>& v)
	{
		vector<double> res = v;
		for(size_t i=0;i<v.size();++i)
		res[i] = v[i]+s;
		return res;
	}
}

void EnclosedMassCalculator::ReCalcEnclosedMass(const double core_mass,const vector<pair<double,double> >& mass_radius_list,
const vector<double>& sample_radii,const double section2shell)
{
	interpolator_=Interpolator(sample_radii,add_all(core_mass,mult_all(section2shell,calc_mass_in_shells(mass_radius_list,sample_radii))));
}

double EnclosedMassCalculator::operator()(double radius) const
{
	return interpolator_(radius);
}

AccelerationCalculator::AccelerationCalculator
(const double gravitation_constant,
const EnclosedMassCalculator& emc):
gravitation_constant_(gravitation_constant), emc_(emc) {}

Vector2D AccelerationCalculator::operator()
(const Tessellation& tess,
const vector<ComputationalCell>& cells,
const vector<Extensive>& /*fluxes*/,
const double /*time*/,
const int point,
TracerStickerNames const& tracerstickernames)const
{
	if(safe_retrieve(cells[static_cast<size_t>(point)].stickers,tracerstickernames.sticker_names,string("ghost")))
      return Vector2D(0,0);
	Vector2D r = tess.GetCellCM(point);
	const double radius = abs(r);
	const double m = emc_(radius);
	return (-1)*gravitation_constant_*m*r/pow(radius,3);
}

CoreAtmosphereGravity::CoreAtmosphereGravity
(const double core_mass,
const vector<double>& sample_radii,
const double gravitation_constant,
const pair<double,double>& sector_angles):
core_mass_(core_mass),
sample_radii_(sample_radii),
gravitation_constant_(gravitation_constant),
section2shell_
(2./(cos(sector_angles.first)-cos(sector_angles.second))),
emc_(EnclosedMassCalculator()),
acc_(AccelerationCalculator(gravitation_constant,emc_))  {}

vector<Extensive> CoreAtmosphereGravity::operator()
(const Tessellation& tess,
const PhysicalGeometry& pg,
const CacheData& cd,
const vector<ComputationalCell>& cells,
const vector<Extensive>& fluxes,
const vector<Vector2D>& point_velocities,
const double time,
const TracerStickerNames& tsn) const
{
	emc_.ReCalcEnclosedMass(core_mass_,calc_mass_radius_list(tess,cells,tsn,cd),sample_radii_,section2shell_);
	ConservativeForce cforce(acc_,true);
	return cforce(tess,pg,cd,cells,fluxes,point_velocities,time,tsn);
}

double CoreAtmosphereGravity::getCoreMass(void) const
{
	return core_mass_;
}

const vector<double>& CoreAtmosphereGravity::getSampleRadii(void) const
{
	return sample_radii_;
}

double CoreAtmosphereGravity::getSection2Shell(void) const
{
	return section2shell_;
}

double CoreAtmosphereGravity::getGravitationConstant(void) const
{
	return gravitation_constant_;
}

