#ifndef LAZY_EXTENSIVE_UPDATER_HPP
#define LAZY_EXTENSIVE_UPDATER_HPP 1

#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"

class LazyExtensiveUpdater: public ExtensiveUpdater
{
public:

  LazyExtensiveUpdater(void);

  void operator()
  (const vector<Extensive>& fluxes,
   const PhysicalGeometry& /*pg*/,
   const Tessellation& tess,
   const double dt,
   const CacheData& cd,
   const vector<ComputationalCell>& cells,
   vector<Extensive>& extensive,
   const TracerStickerNames& tsn) const;
};

#endif // LAZY_EXTENSIVE_UPDATER_HPP
