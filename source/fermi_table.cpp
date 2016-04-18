#include <iostream>
#include "fermi_table.hpp"
#include "source/misc/utils.hpp"

extern "C" {
  void init_tabular_(const char* eos_tab_file);

  void eos_fermi_(int* keyeos,
		  int* im_gas,
		  int* im_photons,
		  int* im_Coulomb,
		  double* rho_in,
		  double* enr_in,
		  double* tmp_in,
		  double* prs_in,
		  double* ent_in,
		  double* abar,
		  double* zbar,
		  double* chem_pot,
		  double* dpdro,
		  double* dpde,
		  double* dedro,
		  double* dedt,
		  double* sounds,
		  int* keyerror);
}

FermiTable::FermiTable
(const string& tab_file,
 const int im_gas,
 const int im_photons,
 const int im_coulomb,
 const boost::container::flat_map<string,pair<double,double> >& atomic_properties):
  im_gas_(im_gas),
  im_photons_(im_photons),
  im_coulomb_(im_coulomb),
  atomic_properties_(atomic_properties)
{
  assert(tab_file.size()<80);
  init_tabular_(tab_file.c_str());
}

double FermiTable::dt2paz(double density, double temperature,
			  std::pair<double,double> aap) const
{
  return calcSingleThermoVar(std::pair<double,double ThermodynamicVariables::*>
			     (density,&ThermodynamicVariables::density),
			     std::pair<double,double ThermodynamicVariables::*>
			     (temperature,&ThermodynamicVariables::temperature),
			     aap,
			     &ThermodynamicVariables::pressure);
}

double FermiTable::dt2p
(double density, 
 double temperature,
 const vector<double>& tracer_values,
 const vector<string>& tracer_names) const
{
  return dt2paz
    (density, 
     temperature, 
     calcAverageAtomicProperties
     (tracer_values,
      tracer_names));
}

double FermiTable::deaz2p(double density, double energy, std::pair<double,double> aap) const
{
  return calcSingleThermoVar(std::pair<double,double ThermodynamicVariables::*>
			     (density,&ThermodynamicVariables::density),
			     std::pair<double,double ThermodynamicVariables::*>
			     (energy,&ThermodynamicVariables::energy),
			     aap,
			     &ThermodynamicVariables::pressure);
}

double FermiTable::dpaz2e(double density, double pressure,
			  std::pair<double, double> aap) const
{
  return calcSingleThermoVar(std::pair<double,double ThermodynamicVariables::*>
			     (density,&ThermodynamicVariables::density),
			     std::pair<double,double ThermodynamicVariables::*>
			     (pressure,&ThermodynamicVariables::pressure),
			     aap,
			     &ThermodynamicVariables::energy);
}

double FermiTable::dpaz2t(double density, double pressure,
			  pair<double,double> aap) const
{
  return calcSingleThermoVar(std::pair<double,double ThermodynamicVariables::*>
			     (density,&ThermodynamicVariables::density),
			     std::pair<double,double ThermodynamicVariables::*>
			     (pressure,&ThermodynamicVariables::pressure),
			     aap,
			     &ThermodynamicVariables::temperature);
}

double FermiTable::deaz2c(double density, double energy,
			  std::pair<double, double> aap) const
{
  return calcSingleThermoVar(std::pair<double,double ThermodynamicVariables::*>
			     (density,&ThermodynamicVariables::density),
			     std::pair<double,double ThermodynamicVariables::*>
			     (energy,&ThermodynamicVariables::energy),
			     aap,
			     &ThermodynamicVariables::sound_speed);
}

double FermiTable::de2c
(double density, 
 double energy, 
 const vector<double>& tracer_values,
 const vector<string>& tracer_names) const
{
  return deaz2c
    (density,
     energy,
     calcAverageAtomicProperties
     (tracer_values,
      tracer_names));
}

double FermiTable::de2p
(double density, 
 double energy, 
 const vector<double>& tracer_values,
 const vector<string>& tracer_names) const
{
  return deaz2p
    (density, 
     energy,
     calcAverageAtomicProperties
     (tracer_values,
      tracer_names));
}

double FermiTable::dpaz2c(double density, double pressure,
			std::pair<double, double> aap) const
{
  return calcSingleThermoVar(std::pair<double,double ThermodynamicVariables::*>
			     (density,&ThermodynamicVariables::density),
			     std::pair<double,double ThermodynamicVariables::*>
			     (pressure,&ThermodynamicVariables::pressure),
			     aap,
			     &ThermodynamicVariables::sound_speed);
}

double FermiTable::dp2c
(double density, 
 double pressure, 
 const vector<double>& tracer_values,
 const vector<string>& tracer_names) const
{
  return dpaz2c
    (density,
     pressure,
     calcAverageAtomicProperties
     (tracer_values,
      tracer_names));
}

double FermiTable::dp2e
(double density,
 double pressure,
 const vector<double>& tracer_values,
 const vector<string>& tracer_names) const
{
  return dpaz2e
    (density,
     pressure,
     calcAverageAtomicProperties
     (tracer_values,
      tracer_names));
}

double FermiTable::dp2t
(double density, 
 double pressure, 
 const vector<double>& tracer_values,
 const vector<string>& tracer_names) const
{
  return dpaz2t
    (density,
     pressure,
     calcAverageAtomicProperties
     (tracer_values,
      tracer_names));
}

double FermiTable::dt2e(double density, double temperature,
			std::pair<double,double> aap) const
{
  return calcSingleThermoVar(std::pair<double,double ThermodynamicVariables::*>
			     (density,&ThermodynamicVariables::density),
			     std::pair<double,double ThermodynamicVariables::*>
			     (temperature,&ThermodynamicVariables::temperature),
			     aap,
			     &ThermodynamicVariables::energy);
}

namespace {
  template<class T> bool compair(const std::pair<T,T>& lhs,
				 const std::pair<T,T>& rhs)
  {
    return ((lhs.first==rhs.first)&&(lhs.second==rhs.second))||
      ((lhs.first==rhs.second)&&(lhs.second==rhs.first));
  }

  typedef std::pair<double FermiTable::ThermodynamicVariables::*,
		    double FermiTable::ThermodynamicVariables::*> tvar_pair;

  FermiTable::Mode determine_mode(const tvar_pair& vars)
  {
    if(compair
       (vars, tvar_pair(&FermiTable::ThermodynamicVariables::density,
			&FermiTable::ThermodynamicVariables::energy)))
      return FermiTable::rho_enr;
    else if(compair
	    (vars, tvar_pair(&FermiTable::ThermodynamicVariables::density,
			     &FermiTable::ThermodynamicVariables::temperature)))
      return FermiTable::rho_tmp;
    else if(compair
	    (vars, tvar_pair(&FermiTable::ThermodynamicVariables::density,
			     &FermiTable::ThermodynamicVariables::pressure)))
      return FermiTable::rho_prs;
    else if(compair
	    (vars, tvar_pair(&FermiTable::ThermodynamicVariables::temperature,
			     &FermiTable::ThermodynamicVariables::pressure)))
      return FermiTable::tmp_prs;
    else if(compair
	    (vars, tvar_pair(&FermiTable::ThermodynamicVariables::energy,
			     &FermiTable::ThermodynamicVariables::pressure)))
      return FermiTable::enr_psr;
    else
      throw "Unsupported combination";
  }		 
}

FermiTable::ThermodynamicVariables::ThermodynamicVariables(void):
  density(1e7), energy(1e7), pressure(1e7), temperature(1e7), entropy(1e7), sound_speed(1e7) {}

double FermiTable::calcSingleThermoVar(std::pair<double,double ThermodynamicVariables::*> input_1,
				       std::pair<double,double ThermodynamicVariables::*> input_2,
				       std::pair<double,double> aap,
				       double ThermodynamicVariables::* output_var) const
{
  ThermodynamicVariables tv;
  tv.*input_1.second = input_1.first;
  tv.*input_2.second = input_2.first;
  calcThermoVars(determine_mode(tvar_pair(input_1.second, input_2.second)),
		 aap, tv);
  return tv.*output_var;
}

void FermiTable::calcThermoVars(Mode mode,
				std::pair<double, double> aap,
				ThermodynamicVariables& tv) const
{
  int keyte = mode;
  double chemical_potential = 0;
  double dpdro = 0;
  double dpde = 0;
  double dedro = 0;
  double dedt = 0;
  int keyerr = 0;
  eos_fermi_(&keyte,
	     &im_gas_,
	     &im_photons_,
	     &im_coulomb_,
	     &tv.density,
	     &tv.energy,
	     &tv.temperature,
	     &tv.pressure,
	     &tv.entropy,
	     &aap.first,
	     &aap.second,
	     &chemical_potential,
	     &dpdro,
	     &dpde,
	     &dedro,
	     &dedt,
	     &tv.sound_speed,
	     &keyerr);
}

std::pair<double,double> FermiTable::calcAverageAtomicProperties
(const vector<double>& tracer_values,
 const vector<string>& tracer_names) const
{
  double total = 0;
  double aa = 0;
  double zz = 0;
  for(boost::container::flat_map<string,std::pair<double,double> >::const_iterator it=
	atomic_properties_.begin();
      it!=atomic_properties_.end();++it){
    //    assert(tracers.count(it->first)==1);
    const double mass_frac = 
      safe_retrieve
      (tracer_values,
       tracer_names,
       it->first);
    const double A = it->second.first;
    const double Z = it->second.second;
    total += mass_frac;
    aa += mass_frac/A;
    zz += mass_frac*Z/A;
  }
  return std::pair<double,double>(total/aa,zz/aa);
}

double FermiTable::dp2s
(double /*density*/, 
 double /*pressure*/,
 const vector<double>& /*tracer_values*/,
 const vector<string>& /*tracer_names*/) const
{
  throw "Method not implemented";
}

double FermiTable::sd2p
(double /*entropy*/,
 double /*density*/,
 const vector<double>& /*tracer_values*/,
 const vector<string>& /*tracer_names*/) const
{
  throw "Method not implemented";
}

const boost::container::flat_map<string,pair<double,double> >& FermiTable::getAtomicProperties(void) const
{
  return atomic_properties_;
}
