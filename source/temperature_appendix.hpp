#ifndef TEMPERATURE_APPENDIX_HPP
#define TEMPERATURE_APPENDIX_HPP 1

#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "fermi_table.hpp"

class TemperatureAppendix: public DiagnosticAppendix
{
public:

  TemperatureAppendix(const FermiTable& eos);

  string getName(void) const;

  vector<double> operator()(const hdsim& sim) const;

private:
  const FermiTable& eos_;
};

#endif // TEMPERATURE_APPENDIX_HPP
