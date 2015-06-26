#include "nuclear_burn.hpp"

extern "C" {
  void burn_step_(int* indexeos,
		  double* density,
		  double* energy,
		  double* tburn,
		  double* xn,
		  double* atomw,
		  double* atomn,
		  double* dedtmp,
		  int* matters,
		  double* dt,
		  double* qrec,
		  int* nse,
		  double* tmp_nse,
		  int* key_done,
		  char* screen_type);
}

namespace {

  pair<double,vector<double> > burn_step_wrapper(double density,
						 double energy,
						 double tburn,
						 vector<double> xn,
						 pair<double,double> az,
						 double dt)
  {
    int indexeos = 0;
    double dedtmp = 0;
    int matters = static_cast<int>(xn.size());
    double qrec = 0;
    int nse = 0;
    double tmp_nse = 1e10;
    char screen_type[80] = "default";
    int key_done = 0;
    burn_step_(&indexeos,
	       &density,
	       &energy,
	       &tburn,
	       &xn[0],
	       &az.first,
	       &az.second,
	       &dedtmp,
	       &matters,
	       &dt,
	       &qrec,
	       &nse,
	       &tmp_nse,
	       &key_done,
	       screen_type);
    assert(key_done==1);
    return pair<double,vector<double> >(qrec,xn);
  }

  template<class S, class T> const T& safe_retrieve(const map<S,T>& m,
						    const S& s)
  {
    /*
    map<S,T>::const_iterator it = m.find(s);
    assert(it!=m.end());
    return it->second;
    */
    assert(m.find(s)!=m.end());
    return m.find(s)->second;
  }

  vector<double> serialize_tracers(const map<string,double>& tracers)
  {
    vector<double> res;
    for(map<string,double>::const_iterator it=tracers.begin();
	it!=tracers.end();
	++it)
      res.push_back(it->second);
    return res;
  }

  map<string,double> reassemble_tracers(const vector<double>& compositions,
					const map<string,double>& old)
  {
    map<string,double> res;
    for(pair<size_t,map<string,double>::const_iterator> index(0,old.begin());
	index.second!=old.end();
	++index.first,++index.second)
      res[index.second->first] = compositions[index.first];
    return res;
  }
}

NuclearBurn::NuclearBurn
(const map<string,pair<double,double> >& atomic_properties,
 const string& ignore_label,
 const FermiTable& eos):
  atomic_properties_(atomic_properties),
  t_prev_(0),
  ignore_label_(ignore_label),
  eos_(eos) {}

void NuclearBurn::operator()(hdsim& sim)
{
  const double dt = sim.getTime() - t_prev_;
  vector<ComputationalCell>& cells = sim.getAllCells();
  for(size_t i=0;i<cells.size();++i){
    ComputationalCell& cell = cells[i];
    if(safe_retrieve(cell.stickers,ignore_label_))
      continue;
    //    const double temperature = eos_.
  }
  sim.recalculateExtensives();
}