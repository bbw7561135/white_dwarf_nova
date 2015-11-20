#ifndef REFLECTIVE_GHOST_THROUGHOUT_HPP
#define REFLECTIVE_GHOST_THROUGHOUT_HPP 1

#include "source/newtonian/two_dimensional/GhostPointGenerator.hpp"

using std::pair;
using std::string;

class ReflectiveGhostThroughout: public GhostPointGenerator
{
public:

  ReflectiveGhostThroughout
  (const string& ghost);

  boost::container::flat_map<size_t, ComputationalCell>
  operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   double time) const;

  pair<ComputationalCell,ComputationalCell>
  GetGhostGradient
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   const vector<pair<ComputationalCell,ComputationalCell> >& gradients,
   size_t ghost_index,
   double time,
   const Edge& edge) const;

private:
  const string ghost_;
};

#endif // REFLECTIVE_GHOST_THROUGHOUT_HPP
