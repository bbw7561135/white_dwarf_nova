/*! \file fermi_table.hpp
  \author Almog Yalinewich
  \brief Equation of state for a degenrate gas
 */

#ifndef FERMI_TABLE_HPP
#define FERMI_TABLE_HPP 1

#include "source/newtonian/common/equation_of_state.hpp"

#include <string>
#include <cassert>
#include <map>
#include <boost/container/flat_map.hpp>

using std::string;
using std::map;
using std::pair;

//! \brief Equation of state of a degenerate gas
class FermiTable: public EquationOfState
{
public:

  /*! \brief Class constructor
    \param tab_file Name of table file
    \param im_gas Gas contribution
    \param im_photons Radiation contribution
    \param im_coulomb Electrostatic contribution
	\param fast_calc Calculates everything from temperature that is given in tracers
   */
  FermiTable
  (const string& tab_file,
   const int im_gas,
   const int im_photons,
   const int im_coulomb,
   const boost::container::flat_map<string,pair<double,double> >& atomic_properties,
   bool fast_calc = true);

  /*! \brief Calculates the pressure
    \param density Density
    \param temperature Temperature
    \param aap First argument is the average weight number, and second is atomic number
    \return Pressure
   */
  double dt2paz
  (double density, 
   double temperature, 
   std::pair<double,double> aap) const;

  double dt2p
  (double density, 
   double temperature,
   const tvector& tracer_values,
   const vector<string>& tracer_names) const;

  /*! \brief Calculates the Energy
    \param density Density
    \param temperature Temperature
    \param aap First argument is the average weight number, and second is atomic number
    \return Energy
   */
  double dt2e(double density, double temperature, std::pair<double,double> aap) const;

  /*! \brief Calculates the temperature
    \param density Density
    \param energy Energy
    \param aap Average atomic properties
    \return Temperature
   */
  double de2t(double density, double energy,  const vector<double>& tracer_values,
   const vector<string>& tracer_names) const;

  /*! \brief Calculates the pressure
    \param density Density
    \param energy Energy
    \param aap Average atomic properties
    \return Pressure
   */
  double deaz2p(double density, double energy, std::pair<double,double> aap) const;

  double de2p
  (double density, 
   double energy, 
   const vector<double>& tracer_values,
   const vector<string>& tracer_names) const;

  /*! \brief Calculates the energy
    \param density Density
    \param pressure Pressure
    \param aap Average atomic properties
    \return Energy
   */
  double dpaz2e(double density, double pressure, std::pair<double, double> aap) const;

  /*! \brief Calculates the speed of sound
    \param density Density
    \param energy Energy
    \param aap Average atomic properties
    \return Speed of sound
   */
  double deaz2c(double density, double energy, std::pair<double, double> aap) const;

  double de2c
  (double density, 
   double energy, 
   const vector<double>& tracer_values,
   const vector<string>& tracer_name) const;

  /*! \brief Calculates the speed of sound
    \param density Density
    \param pressure Pressure
    \param aap Average atomic properties
    \return Speed of sound
   */
  double dpaz2c(double density, double pressure, std::pair<double, double> aap) const;

  double dpaz2t(double density, double pressure, pair<double,double> aap) const;
  
  double dtaz2c(double density, double temperature,pair<double, double> aap) const;

  double dp2c
  (double density, 
   double pressure, 
   const vector<double>& tracer_values,
   const vector<string>& tracer_names) const;
   
   double dt2c
  (double density, 
   double temperature, 
   const vector<double>& tracer_values,
   const vector<string>& tracer_names) const;

  double dp2e
  (double density, 
   double pressure, 
   const vector<double>& tracer_values,
   const vector<string>& tracer_names) const;

  double dp2t
  (double density,
   double pressure,
   const vector<double>& tracer_values,
   const vector<string>& tracer_names) const;

  double dp2s
  (double density, 
   double pressure, 
   const vector<double>& tracer_values,
   const vector<string>& tracer_names) const;

  double sd2p
  (double entropy, 
   double density, 
   const vector<double>& tracer_values,
   const vector<string>& tracer_names) const;

  enum Mode{
    rho_enr,
    rho_tmp,
    rho_prs,
    tmp_prs,
    enr_psr
  };

  class ThermodynamicVariables
  {
  public:

    double density;
    double energy;
    double pressure;
    double temperature;
    double entropy;
    double sound_speed;

    ThermodynamicVariables(void);
  };

  double calcSingleThermoVar(std::pair<double,double ThermodynamicVariables::*> input_1,
			     std::pair<double,double ThermodynamicVariables::*> input_2,
			     std::pair<double,double> aap,
			     double ThermodynamicVariables::* output_var) const;
				 

  void calcThermoVars(Mode mode,
		      std::pair<double, double> aap,
		      ThermodynamicVariables& tv) const;

  std::pair<double,double> calcAverageAtomicProperties
  (const vector<double>& tracer_values,
   const vector<string>& tracer_names) const;

  const boost::container::flat_map<string,pair<double,double> >& getAtomicProperties(void) const;

  void index_check(vector<string> const& names) const;
private:
  mutable int im_gas_;
  mutable int im_photons_;
  mutable int im_coulomb_;
  const boost::container::flat_map<string,std::pair<double,double> > atomic_properties_;
  const bool fast_calc_;
  mutable int temperature_index_;
  void index_check(vector<string> const& names);
};

#endif // FERMI_TABLE_HPP
