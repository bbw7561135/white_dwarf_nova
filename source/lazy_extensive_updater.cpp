#include "lazy_extensive_updater.hpp"
#include "bracketed.hpp"
#include <string>
#include "source/misc/utils.hpp"

using std::string;

LazyExtensiveUpdater::LazyExtensiveUpdater(void) {}

void LazyExtensiveUpdater::operator()
(const vector<Extensive>& fluxes,
 const PhysicalGeometry& /*pg*/,
 const Tessellation& tess,
 const double dt,
 const CacheData& cd,
 const vector<ComputationalCell>& cells,
 vector<Extensive>& extensive) const
{
  const vector<Edge>& edge_list = tess.getAllEdges();
  for(size_t i=0;i<edge_list.size();++i){
    const Edge& edge = edge_list.at(i);
    const Extensive delta = dt*cd.areas[i]*fluxes.at(i);
    if(bracketed(0,edge.neighbors.first,tess.GetPointNo()) &&
       !safe_retrieve
       (cells.at(static_cast<size_t>(edge.neighbors.first)).stickers,string("ghost")))
      extensive.at(static_cast<size_t>(edge.neighbors.first)) -= delta;
    if(bracketed(0,edge.neighbors.second,tess.GetPointNo()) &&
       !safe_retrieve
       (cells.at(static_cast<size_t>(edge.neighbors.second)).stickers,string("ghost")))
      extensive.at(static_cast<size_t>(edge.neighbors.second)) += delta;
  }
}
