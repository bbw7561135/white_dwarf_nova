#ifndef NUCLEAR_BURN_HPP
#define NUCLEAR_BURN_HPP 1

#include <map>
#include <string>
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "fermi_table.hpp"

using std::map;
using std::string;
using std::pair;

class NuclearBurn: public Manipulate
{
public:
  NuclearBurn(const string& rfile,
	      const string& ignore_label,
	      const FermiTable& eos,
	      const string& ehf);

  void operator()(hdsim& sim);

  double suggestInverseTimeStep(void) const
  {
    return idt_;
  }

  ~NuclearBurn(void);
  
private:

  mutable double idt_;
  mutable double t_prev_;
  const string ignore_label_;
  const FermiTable& eos_;
  const vector<string> isotope_list_;
  const string energy_history_fname_;
  mutable vector<pair<double, double> > energy_history_;
};

#endif // NUCLEAR_BURN_HPP
