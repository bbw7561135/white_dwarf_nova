#include "nuclear_burn.hpp"
#include "source/misc/vector_initialiser.hpp"
#include <fstream>

extern "C" {

  void initnet_(const char* rfile);

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
    if(key_done!=1){
      std::ofstream f("burn_step_error_report.txt");
      f << "density = " << density << "\n";
      f << "energy = " << energy << "\n";
      f << "temperature = " << tburn << "\n";
      f << "atomic weight = " << az.first << "\n";
      f << "atomic number = " << az.second << "\n";
      f << "dt = " << dt << "\n";
      f.close();
      assert(key_done==1);
    }
    return pair<double,vector<double> >(qrec,xn);
  }

  vector<double> serialize_tracers
  (const boost::container::flat_map<string,double>& tracers,
   const vector<string>& isotope_list)
  {
    vector<double> res;
    for(size_t i=0;i<isotope_list.size();++i){
      res.push_back(safe_retrieve(tracers,isotope_list.at(i)));
    }
    return res;
  }

  boost::container::flat_map<string,double> reassemble_tracers
  (const vector<double>& compositions,
   const vector<string>& isotope_list)
  {
    assert(compositions.size()==isotope_list.size());
    boost::container::flat_map<string,double> res;
    for(size_t i=0;i<compositions.size();++i)
      res[isotope_list[i]] = compositions[i];
    return res;
  }
}

NuclearBurn::NuclearBurn
(const string& rfile,
 const string& ignore_label,
 const FermiTable& eos,
 const string& ehf):
  t_prev_(0),
  ignore_label_(ignore_label),
  eos_(eos),
  isotope_list_(VectorInitialiser<string>("He4")
		("C12")
		("O16")
		("Ne20")
		("Mg24")
		("Si28")
		("S32")
		("Ar36")
		("Ca40")
		("Ti44")
		("Cr48")
		("Fe52")
		("Ni56")()),
  energy_history_fname_(ehf),
  energy_history_()
{
  initnet_(rfile.c_str());
}

void NuclearBurn::operator()(hdsim& sim)
{
  const double dt = sim.getTime() - t_prev_;
  t_prev_ = sim.getTime();
  double total = 0;
  vector<ComputationalCell>& cells = sim.getAllCells();
  for(size_t i=0;i<cells.size();++i){
    ComputationalCell& cell = cells[i];
    if(safe_retrieve(cell.stickers,ignore_label_))
      continue;
    const double temperature = eos_.dp2t(cell.density,
					 cell.pressure,
					 cell.tracers);
    const double energy = eos_.dp2e(cell.density,
				    cell.pressure,
				    cell.tracers);
    const pair<double,vector<double> > qrec_tracers =
      burn_step_wrapper(cell.density,energy,temperature,
			serialize_tracers(cell.tracers,
					  isotope_list_),
			eos_.calcAverageAtomicProperties(cell.tracers),
			dt);
    total += dt*qrec_tracers.first;
    const double new_energy = energy + dt*qrec_tracers.first;
    cell.tracers = reassemble_tracers(qrec_tracers.second,isotope_list_);
    cell.pressure = eos_.de2p(cell.density, new_energy, cell.tracers);
  }
  sim.recalculateExtensives();
  energy_history_.push_back
    (pair<double,double>(sim.getTime(),total));
}

NuclearBurn::~NuclearBurn(void)
{
  std::ofstream f(energy_history_fname_.c_str());
  for(size_t i=0;i<energy_history_.size();++i)
    f << energy_history_[i].first << " "
      << energy_history_[i].second << std::endl;
  f.close();
}
