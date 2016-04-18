#include "create_pressure_reference.hpp"

vector<double> create_pressure_reference
(const FermiTable& eos,
 const InitialData& id)
{
  vector<double> res(id.radius_list.size());
  vector<string> tracer_names;
  for(boost::container::flat_map<string,vector<double> >::const_iterator it=
	id.tracers_list.begin();
      it!=id.tracers_list.end();
      ++it)
    tracer_names.push_back(it->first);
  for(size_t i=0;i<id.density_list.size();++i){
    vector<double> tracer_values;;
    for(boost::container::flat_map<string,vector<double> >::const_iterator it=
	  id.tracers_list.begin();
	it!=id.tracers_list.end();
	++it)
      tracer_values.push_back((it->second).at(i));
    res[i] = eos.dt2p(id.density_list.at(i),
		      id.temperature_list.at(i),
		      tracer_values,
		      tracer_names);
  }
  return res;
}
