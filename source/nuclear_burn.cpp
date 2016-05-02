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

  pair<double,vector<double> > burn_step_wrapper
  (double density,
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

namespace
{
	vector<double> trim_tracers(vector<double> const& tracers,vector<size_t> const& indeces)
	{
		size_t N=indeces.size();
		vector<double> res(N);
		for(size_t i=0;i<N;++i)
			res[i]=tracers[indeces[i]];
		return res;
	}
	
	void update_tracers(vector<double> &tracers,vector<size_t> const& indeces,vector<double> const& trimed_vector)
	{
		size_t N=indeces.size();
		for(size_t i=0;i<N;++i)
			tracers[indeces[i]]=trimed_vector[i];
	}
}

void NuclearBurn::operator()(hdsim& sim)
{
  const double dt = sim.getTime() - t_prev_;
  t_prev_ = sim.getTime();
  double total = 0;
  vector<ComputationalCell>& cells = sim.getAllCells();
  vector<Extensive> extensives = sim.getAllExtensives();
  size_t Nloop = static_cast<size_t>(sim.getTessellation().GetPointNo());
  
  vector<size_t> fortran_indeces(isotope_list_.size());
  TracerStickerNames tracerstickernames =  sim.GetTracerStickerNames();
  for(size_t i=0;i<fortran_indeces.size();++i)
  {
	  vector<string>::const_iterator it = binary_find(tracerstickernames.tracer_names.begin(),
			tracerstickernames.tracer_names.end(), isotope_list_[i]);
	  assert(it != tracerstickernames.tracer_names.end());
	  fortran_indeces.at(i) = static_cast<size_t>(it - tracerstickernames.tracer_names.begin());
  }
  vector<string>::const_iterator it2 = binary_find(tracerstickernames.tracer_names.begin(),
			tracerstickernames.tracer_names.end(), string("Temperature"));
  assert(it2 != tracerstickernames.tracer_names.end());
  size_t temperature_index = static_cast<size_t>(it2-tracerstickernames.tracer_names.begin());
  
  for(size_t i=0;i<Nloop;++i)
  {
    ComputationalCell& cell = cells[i];
    if(safe_retrieve
       (cell.stickers,
	sim.GetTracerStickerNames().sticker_names,
	ignore_label_))
      continue;
    const double energy = (extensives[i].energy-0.5*ScalarProd(extensives[i].momentum,extensives[i].momentum)/extensives[i].mass)
		/extensives[i].mass;
	 	   
    const pair<double,vector<double> > qrec_tracers =
      burn_step_wrapper
      (cell.density,energy,cell.tracers[temperature_index],
		trim_tracers(cell.tracers,fortran_indeces),
       //cell.tracers,
       eos_.calcAverageAtomicProperties
       (cell.tracers,
	sim.GetTracerStickerNames().tracer_names),
       dt);
    total += dt*qrec_tracers.first;
    const double new_energy = energy + dt*qrec_tracers.first;
	
    update_tracers(cell.tracers,fortran_indeces,qrec_tracers.second);
	const double temperature_new = eos_.de2t(cell.density, new_energy, cell.tracers, sim.GetTracerStickerNames().tracer_names);
 	cell.tracers[temperature_index]=temperature_new;
    cell.pressure = eos_.de2p
      (cell.density,
       new_energy,
       cell.tracers,
       sim.GetTracerStickerNames().tracer_names);
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
