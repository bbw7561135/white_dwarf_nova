#include <cmath>
#include "calc_init_cond.hpp"
#include "create_pressure_reference.hpp"
#include "vector_io.hpp"
#include "interpolator.hpp"
#include "source/misc/utils.hpp"

pair<TracerStickerNames, vector<ComputationalCell> >
calc_init_cond(const Tessellation& tess,
	       const FermiTable& eos,
	       const InitialData& id,
	       const Shape2D& cd)
{
  save_txt("pressure_reference.txt",create_pressure_reference(eos,id));
  vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
  const Interpolator density_interpolator(id.radius_mid,
					  id.density_list);
  const Interpolator temperature_interpolator(id.radius_mid,
					      id.temperature_list);
  const Interpolator velocity_interpolator(id.radius_list,
					   id.velocity_list);
  boost::container::flat_map<string,Interpolator*> tracer_intepolators;
  for(boost::container::flat_map<string,vector<double> >::const_iterator it=
	id.tracers_list.begin();
      it!=id.tracers_list.end(); ++it)
    tracer_intepolators[it->first] = new Interpolator(id.radius_mid,
						      it->second);
  TracerStickerNames res_2;
  res_2.sticker_names.push_back("ghost");
  for(boost::container::flat_map<string,Interpolator*>::const_iterator it=
	tracer_intepolators.begin();
      it!=tracer_intepolators.end();
      ++it)
    res_2.tracer_names.push_back(it->first);
  sort(res_2.tracer_names.begin(), 
       res_2.tracer_names.end());
  for(size_t i=0;i<res.size();++i){
    res.at(i).density = id.density_list.back();
    res.at(i).velocity = Vector2D(0,0);
    res.at(i).stickers[0] = true;
    res.at(i).tracers.resize(res_2.tracer_names.size());
    for(boost::container::flat_map<string,Interpolator*>::const_iterator it=
	  tracer_intepolators.begin();
	it!=tracer_intepolators.end();
	++it)
      safe_retrieve
	(res.at(i).tracers,
	 res_2.tracer_names,
	 it->first) = 0;
    safe_retrieve
      (res.at(i).tracers,
       res_2.tracer_names,
       string("He4")) = 1;
    res.at(i).pressure = eos.dt2p(res.at(i).density,
				  id.temperature_list.back(),
				  res.at(i).tracers,
				  res_2.tracer_names);
    const Vector2D r = tess.GetCellCM(static_cast<int>(i));
    const double radius = abs(r);
    if(!cd(r))
      continue;
    safe_retrieve
      (res.at(i).stickers,
       res_2.sticker_names,
       string("ghost")) = false;
    const double density = density_interpolator(radius);
    const double temperature = temperature_interpolator(radius);
    const double velocity = velocity_interpolator(radius);
    for(boost::container::flat_map<string,Interpolator*>::const_iterator it=
	  tracer_intepolators.begin();
	it!=tracer_intepolators.end();
	++it)
      safe_retrieve
	(res.at(i).tracers,
	 res_2.tracer_names,
	 it->first) = (*(it->second))(radius);
    const double pressure = eos.dt2p
      (density, 
       temperature, 
       res.at(i).tracers,
       res_2.tracer_names);
    res.at(i).density = density;
    res.at(i).pressure = pressure;
    res.at(i).velocity = r*velocity/radius;
  }
  for(boost::container::flat_map<string,Interpolator*>::iterator it=
	tracer_intepolators.begin();
      it!=tracer_intepolators.end();
      ++it)
    delete it->second;
  return pair<TracerStickerNames, vector<ComputationalCell> >
    (res_2, res);
}
