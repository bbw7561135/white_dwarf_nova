#include "monopole_self_gravity.hpp"
#include "interpolator.hpp"

namespace {
  vector<pair<double,double> > calc_mass_radius_list
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   const CacheData& cd,
   const TracerStickerNames& tsn)
  {
    vector<pair<double, double> > res;
    for(size_t i=0;i<cells.size();++i){
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

  vector<double> calc_mass_in_shells(const vector<pair<double,double> >& mass_radius_list,
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
	  //	  break;
	}
      }
    }
    return res;
  }
}
  

MonopoleSelfGravity::MonopoleSelfGravity
(const vector<double>& sample_radii,
 double gravitation_constant,
 const pair<double,double>& section_angles):
  sample_radii_(sample_radii),
  gravitation_constant_(gravitation_constant),
  section2shell_(2./(cos(section_angles.first)-cos(section_angles.second))) {}

vector<Extensive> MonopoleSelfGravity::operator()
  (const Tessellation& tess,
   const PhysicalGeometry& /*pg*/,
   const CacheData& cd,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& /*fluxes*/,
   const vector<Vector2D>& /*point_velocities*/,
   const double /*t*/,
   const TracerStickerNames& tsn) const
{
  const vector<double> mass_sample =
    calc_mass_in_shells
    (calc_mass_radius_list
     (tess,
      cells,
      cd,
      tsn),
     sample_radii_);
  const Interpolator radius_mass_interp(sample_radii_,
					mass_sample);

  vector<Extensive> res(static_cast<size_t>(tess.GetPointNo()));
  for(size_t i=0;i<res.size();++i){
    if(safe_retrieve
       (cells.at(i).stickers,
	tsn.sticker_names,
	string("ghost")))
      continue;
    const Vector2D r = tess.GetCellCM(static_cast<int>(i));
    const double radius = abs(r);
    const double mass = radius_mass_interp(radius)*section2shell_;
    const Vector2D acc = (-1)*gravitation_constant_*r*mass/pow(radius,3);
    const double volume = cd.volumes[i];
    res[i].mass = 0;
    res[i].momentum = volume*cells[i].density*acc;
    res[i].energy = volume*cells[i].density*ScalarProd(acc,cells[i].velocity);
  }
  return res;
}
