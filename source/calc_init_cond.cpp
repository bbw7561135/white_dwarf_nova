#include <cmath>
#include "calc_init_cond.hpp"
#include "create_pressure_reference.hpp"
#include "vector_io.hpp"
#include "interpolator.hpp"

vector<ComputationalCell> calc_init_cond(const Tessellation& tess,
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
  for(size_t i=0;i<res.size();++i){
    res.at(i).density = id.density_list.back();
    res.at(i).velocity = Vector2D(0,0);
    res.at(i).stickers["ghost"] = true;
    for(boost::container::flat_map<string,Interpolator*>::const_iterator it=
	  tracer_intepolators.begin();
	it!=tracer_intepolators.end();
	++it)
      res.at(i).tracers[it->first] = 0;
    res.at(i).tracers["He4"] = 1;
    res.at(i).pressure = eos.dt2p(res.at(i).density,
				  id.temperature_list.back(),
				  res.at(i).tracers);
    const Vector2D r = tess.GetCellCM(static_cast<int>(i));
    const double radius = abs(r);
    if(!cd(r))
      continue;
    res.at(i).stickers["ghost"] = false;
    const double density = density_interpolator(radius);
    const double temperature = temperature_interpolator(radius);
    const double velocity = velocity_interpolator(radius);
    for(boost::container::flat_map<string,Interpolator*>::const_iterator it=
	  tracer_intepolators.begin();
	it!=tracer_intepolators.end();
	++it)
      res.at(i).tracers[it->first] = (*(it->second))(radius);
    const double pressure = eos.dt2p(density, temperature, res.at(i).tracers);
    res.at(i).density = density;
    res.at(i).pressure = pressure;
    res.at(i).velocity = r*velocity/radius;
  }
  for(boost::container::flat_map<string,Interpolator*>::iterator it=
	tracer_intepolators.begin();
      it!=tracer_intepolators.end();
      ++it)
    delete it->second;
  return res;
}
