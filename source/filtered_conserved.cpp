#include <fstream>
#include "filtered_conserved.hpp"

using namespace std;

FilteredConserved::FilteredConserved(const string& fname):
  data_(), fname_(fname) {}

void FilteredConserved::operator()(const hdsim& sim)
{
  const vector<ComputationalCell>& cells = sim.getAllCells();
  const vector<Extensive>& extensives = sim.getAllExtensives();
  Extensive buf;
  buf.mass = 0;
  buf.momentum = Vector2D(0,0);
  buf.energy = 0;
  for(size_t i=0;i<extensives.size();++i){
    if(safe_retrieve(cells[i].stickers,string("ghost")))
      continue;
    buf += extensives[i];
  }
  data_.push_back
    (pair<double,Extensive>
     (sim.getTime(),buf));
}

FilteredConserved::~FilteredConserved(void)
{
  ofstream f(fname_.c_str());
  for(size_t i=0;i<data_.size();++i)
    f << data_[i].first << " "
      << data_[i].second.mass << " "
      << data_[i].second.momentum.x << " "
      << data_[i].second.momentum.y << " "
      << data_[i].second.energy << endl;
  f.close();
}
