#include "temperature_appendix.hpp"

TemperatureAppendix::TemperatureAppendix(const FermiTable& eos):
  eos_(eos) {}

string TemperatureAppendix::getName(void) const
{
  return "temperature";
}

vector<double> TemperatureAppendix::operator()(const hdsim& sim) const
{
  const vector<ComputationalCell>& cells = sim.getAllCells();
  vector<double> temperatures(cells.size(),1);
  for(size_t i=0;i<cells.size();++i){
    const ComputationalCell& cell = cells[i];
    if(cell.stickers.front())
      continue;
    temperatures[i] = eos_.dp2t
      (cell.density,
       cell.pressure,
       cell.tracers,
       sim.GetTracerStickerNames().tracer_names);
  }
  return temperatures;
}
