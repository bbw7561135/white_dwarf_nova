#include "reflective_ghost_throughout.hpp"
#include "source/misc/utils.hpp"

ReflectiveGhostThroughout::ReflectiveGhostThroughout
(const string& ghost): ghost_(ghost) {}

namespace {
  ComputationalCell reverse_normal_velocity
  (const ComputationalCell& c,
   const Edge& e)
  {
    const Vector2D p =
      normalize(e.vertices.second - e.vertices.first);
    ComputationalCell res = c;
    res.velocity = 2*ScalarProd(res.velocity,p)*p - res.velocity;
    return res;
  }

  vector<pair<size_t, size_t> > get_special_interface
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   const string& ghost)
  {
    vector<pair<size_t, size_t> > res;
    const vector<Edge>& edge_list = tess.getAllEdges();
    for(size_t i=0;i<edge_list.size();++i){
      const Edge& edge = edge_list[i];
      if(edge.neighbors.first<0 ||
	 edge.neighbors.second<0 ||
	 edge.neighbors.first>tess.GetPointNo() ||
	 edge.neighbors.second>tess.GetPointNo())
	continue;
      if(safe_retrieve
	 (cells.at
	  (static_cast<size_t>
	   (edge.neighbors.first)).stickers,
	  ghost) &&
	 !safe_retrieve
	 (cells.at
	  (static_cast<size_t>
	   (edge.neighbors.second)).stickers,
	  ghost))
	res.push_back(pair<size_t,size_t>(i,1));
      else if(safe_retrieve
	      (cells.at
	       (static_cast<size_t>
		(edge.neighbors.second)).stickers,
	       ghost) &&
	      !safe_retrieve
	      (cells.at
	       (static_cast<size_t>
		(edge.neighbors.first)).stickers,
	       ghost))
	res.push_back(pair<size_t,size_t>(i,2));
    }
    return res;
  }
}

boost::container::flat_map<size_t, ComputationalCell>
ReflectiveGhostThroughout::operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   double /*time*/) const
{
  boost::container::flat_map<size_t, ComputationalCell> res;
  const vector<pair<size_t, size_t> > ghosts =
    join
    (GetOuterEdgesIndeces(tess),
     get_special_interface
     (tess,
      cells,
      ghost_));
  for(size_t i=0;i<ghosts.size();++i){
    const Edge& edge = tess.GetEdge
      (static_cast<int>(ghosts[i].first));
    const ComputationalCell c =
      reverse_normal_velocity
      (cells.at
       (static_cast<size_t>
	(ghosts[i].second==2 ?
	 edge.neighbors.first :
	 edge.neighbors.second)),
       edge);
    res.insert
      (pair<size_t,ComputationalCell>
       (static_cast<size_t>
	(ghosts[i].second == 1 ?
	 edge.neighbors.first :
	 edge.neighbors.second),
	c));
  }
  return res;
}

pair<ComputationalCell,ComputationalCell>
ReflectiveGhostThroughout::GetGhostGradient
(const Tessellation& tess,
 const vector<ComputationalCell>& /*cells*/,
 const vector<pair<ComputationalCell,ComputationalCell> >& gradients,
 size_t ghost_index,
 double /*time*/,
 const Edge& /*edge*/) const
{
  pair<ComputationalCell,ComputationalCell> grad = 
    gradients
    [static_cast<size_t>
     (tess.GetOriginalIndex
      (static_cast<int>(ghost_index)))];
  const Vector2D normal = 
    normalize
    (tess.GetMeshPoint(static_cast<int>(ghost_index)) - 
     tess.GetMeshPoint
     (tess.GetOriginalIndex(static_cast<int>(ghost_index))));
  grad.first.velocity -= 2*ScalarProd(grad.first.velocity,normal)*normal;
  grad.second.velocity -= 2*ScalarProd(grad.second.velocity, normal)*normal;
  return grad;
}
