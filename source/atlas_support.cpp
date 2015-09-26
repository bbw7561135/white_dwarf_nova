#include "atlas_support.hpp"

AtlasSupport::AtlasSupport(void) {}

namespace {

  class Bracket
  {
  public:

    Bracket(int low, int high):
      low_high_(low,high) {}

    bool operator()(int arg) const
    {
      return arg>=low_high_.first && arg<=low_high_.second;
    }

  private:
    const pair<int,int> low_high_;
  };

  double calc_radius_sqr(const Vector2D& p)
  {
    return ScalarProd(p,p);
  }

  bool point_above_edge(const Vector2D& p,
			const Edge& e)
  {
    const double r2 = calc_radius_sqr(p);
    return r2>calc_radius_sqr(e.vertices.first) &&
      r2>calc_radius_sqr(e.vertices.second);
  }

  vector<size_t> get_lower_layer_indices
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells)
  {
    vector<size_t> res;
    const vector<Edge>& edge_list = tess.getAllEdges();
    const Bracket b(0,tess.GetPointNo());
    for(size_t i=0;i<edge_list.size();++i){
      const Edge& edge = edge_list[i];
      if(!(b(edge.neighbors.first)&&
	   b(edge.neighbors.second)))
	continue;
      const size_t left_index = 
	static_cast<size_t>(edge.neighbors.first);
      const size_t right_index =
	static_cast<size_t>(edge.neighbors.second);
      if(safe_retrieve(cells.at(left_index).stickers,
		       string("ghost")) &&
	 !safe_retrieve(cells.at(right_index).stickers,
			string("ghost")) &&
	 point_above_edge(tess.GetCellCM(edge.neighbors.second),
			  edge))
	res.push_back(right_index);
      else if
	(safe_retrieve(cells.at(right_index).stickers,
		       string("ghost")) &&
	 !safe_retrieve(cells.at(left_index).stickers,
			string("ghost")) &&
	 point_above_edge(tess.GetCellCM(edge.neighbors.first),
			  edge))
	res.push_back(left_index);
    }
    return res;
  }
}

void AtlasSupport::operator()(hdsim& sim)
{
  const Tessellation& tess = sim.getTessellation();
  const vector<ComputationalCell>& cell_list = sim.getAllCells();
  const vector<size_t> index_list = 
    get_lower_layer_indices(tess,cell_list);
  vector<Extensive>& extensive_list = sim.getAllExtensives();
  for(size_t i=0;i<index_list.size();++i){
    const Vector2D r = tess.GetMeshPoint
      (static_cast<int>(index_list[i]));
    extensive_list[index_list[i]].momentum = 
      0*abs(extensive_list[index_list[i]].momentum)*
      r/abs(r);
  }
  sim.recalculatePrimitives();
}
