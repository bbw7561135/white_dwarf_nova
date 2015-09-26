#include "calc_bottom_area.hpp"

namespace {

  class Bracket
  {
  public:

    Bracket(int low,
	    int high):
      low_(low), high_(high) {}

    bool operator()(int n) const
    {
      return n>=low_ && n<=high_;
    }

  private:
    int low_;
    int high_;
  };

  double calc_radius_sqr(const Vector2D& p)
  {
    return ScalarProd(p,p);
  }

  bool point_above_edge(const Vector2D& p,
			const Edge& edge)
  {
    return (calc_radius_sqr(p)>calc_radius_sqr(edge.vertices.first) &&
	    calc_radius_sqr(p)>calc_radius_sqr(edge.vertices.second));
  }

  bool is_bottom_edge(const Edge& edge,
		      const Tessellation& tess,
		      const Bracket& b,
		      const Shape2D& shape)
  {
    if(!b(edge.neighbors.first) || !b(edge.neighbors.second))
      return false;
    if(!shape(tess.GetMeshPoint(edge.neighbors.first))){
      if(!shape(tess.GetMeshPoint(edge.neighbors.second)))
	return false;
      const Vector2D r = tess.GetCellCM(edge.neighbors.second);
      if(point_above_edge(r,edge))
	return true;
      return false;
    }
    if(!shape(tess.GetMeshPoint(edge.neighbors.second))){
      const Vector2D r = tess.GetCellCM(edge.neighbors.first);
      if(point_above_edge(r,edge))
	return true;
      return false;
    }
    return false;
  }
}

double calc_bottom_area(const Tessellation& tess,
			const Shape2D& shape,
			const PhysicalGeometry& pg)
{
  const vector<Edge>& edge_list = tess.getAllEdges();
  const Bracket b(0,tess.GetPointNo());
  double res = 0;
  for(size_t i=0;i<edge_list.size();++i){
    const Edge& edge = edge_list[i];
    if(is_bottom_edge(edge,tess,b,shape))
      res += pg.calcArea(edge);
  }
  return res;
}
