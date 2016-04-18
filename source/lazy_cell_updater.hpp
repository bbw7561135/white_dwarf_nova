#ifndef LAZY_CELL_UPDATER_HPP
#define LAZY_CELL_UPDATER_HPP 1

#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"

class LazyCellUpdater: public CellUpdater
{
public:

  LazyCellUpdater(void);

  vector<ComputationalCell> operator()
  (const Tessellation& /*tess*/,
   const PhysicalGeometry& /*pg*/,
   const EquationOfState& eos,
   vector<Extensive>& extensives,
   const vector<ComputationalCell>& old,
   const CacheData& cd,
   const TracerStickerNames& tsn) const;
};

#endif // LAZY_CELL_UPDATER_HPP
