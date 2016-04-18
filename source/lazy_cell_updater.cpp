#include "lazy_cell_updater.hpp"

LazyCellUpdater::LazyCellUpdater(void) {}

vector<ComputationalCell> LazyCellUpdater::operator()
  (const Tessellation& /*tess*/,
   const PhysicalGeometry& /*pg*/,
   const EquationOfState& eos,
   vector<Extensive>& extensives,
   const vector<ComputationalCell>& old,
   const CacheData& cd,
   const TracerStickerNames& tsn) const
{
  vector<ComputationalCell> res = old;
  for(size_t i=0;i<extensives.size();++i){
    if(safe_retrieve
       (old.at(i).stickers,
	tsn.tracer_names,
	string("ghost")))
      continue;
    const double volume = cd.volumes[i];
    res.at(i).density = extensives.at(i).mass/volume;
    res.at(i).velocity = extensives.at(i).momentum/extensives.at(i).mass;
    const double total_energy = extensives.at(i).energy/extensives.at(i).mass;
    const double kinetic_energy = 0.5*ScalarProd(res.at(i).velocity, res.at(i).velocity);
    const double thermal_energy = total_energy - kinetic_energy;
    res.at(i).tracers = vector<double>
      (extensives.at(i).tracers.size(),0);
    for(size_t j=0;i<extensives.at(j).tracers.size();++j)
      res.at(i).tracers.at(j) = 
	extensives.at(i).tracers.at(j)/
	extensives.at(i).mass;
    res.at(i).pressure = 
      eos.de2p
      (res.at(i).density,
       thermal_energy, 
       res.at(i).tracers,
       tsn.tracer_names);
  }
  return res;
}
