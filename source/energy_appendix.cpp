#include "energy_appendix.hpp"

EnergyAppendix::EnergyAppendix(const FermiTable& eos):
  eos_(eos) {}

string EnergyAppendix::getName(void) const
{
  return "energy";
}

vector<double> EnergyAppendix::operator()(const hdsim& sim) const
{
  const vector<ComputationalCell>& cells = sim.getAllCells();
  vector<double> res(cells.size(), 1);
  for(size_t i=0;i<res.size();++i){
    const ComputationalCell& cell = cells[i];
    if(safe_retrieve(cell.stickers,string("ghost")))
      continue;
    res[i] = eos_.dp2e(cell.density,
		       cell.pressure,
		       cell.tracers);
  }
  return res;
}
