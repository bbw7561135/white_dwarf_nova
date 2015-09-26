#include "inner_bc.hpp"
#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"

namespace {

  Extensive conserved_to_extensive
  (const Conserved& c, const ComputationalCell& cell)
  {
    Extensive res;
    res.mass = c.Mass;
    res.momentum = c.Momentum;
    res.energy = c.Energy;
    for(boost::container::flat_map<string,double>::const_iterator it=
	  cell.tracers.begin();
	it!=cell.tracers.end();++it)
      res.tracers[it->first] = (it->second)*c.Mass;
    return res;
  }
}

InnerBC::InnerBC
(const RiemannSolver& rs,
 const string& ghost,
 const EquationOfState& eos):
  rs_(rs),
  ghost_(ghost),
  eos_(eos) {}

vector<pair<size_t,Extensive> > InnerBC::operator()
  (const Tessellation& tess, 
   const vector<ComputationalCell>& cells) const
{
  vector<pair<size_t,Extensive> > res;
  const vector<Edge>& edge_list = tess.getAllEdges();
  for(size_t i=0;i<edge_list.size();++i){
    const Edge& edge = edge_list[i];
    if(edge.neighbors.first<0 ||
       edge.neighbors.second<0 ||
       edge.neighbors.first>tess.GetPointNo() ||
       edge.neighbors.second>tess.GetPointNo())
      continue;
    if(safe_retrieve(cells.at(static_cast<size_t>(edge.neighbors.first)).stickers,ghost_)){
      if(safe_retrieve(cells.at(static_cast<size_t>(edge.neighbors.second)).stickers,ghost_)){
	res.push_back(pair<size_t,Extensive>(i,Extensive()));
	continue;
      }
      const Vector2D p = 
	normalize(edge.vertices.second - edge.vertices.first);
      const Vector2D n =
	normalize(tess.GetMeshPoint(edge.neighbors.second)-
		  tess.GetMeshPoint(edge.neighbors.first));
      const ComputationalCell& right_cell = 
	cells.at(static_cast<size_t>(edge.neighbors.second));
      const Primitive right = 
	convert_to_primitive(right_cell,eos_);
      const Conserved c = rotate_solve_rotate_back
	(rs_,
	 reflect(right,p),
	 right,
	 0,n,p);
      res.push_back
	(pair<size_t,Extensive>
	 (i,conserved_to_extensive
	  (c,right_cell)));
    }
    if(safe_retrieve(cells.at(static_cast<size_t>(edge.neighbors.second)).stickers,ghost_)){
      const Vector2D p = 
	normalize(edge.vertices.second - edge.vertices.first);
      const Vector2D n =
	normalize(tess.GetMeshPoint(edge.neighbors.second)-
		  tess.GetMeshPoint(edge.neighbors.first));
      const ComputationalCell& left_cell = 
	cells.at(static_cast<size_t>(edge.neighbors.first));
      const Primitive left = 
	convert_to_primitive(left_cell,eos_);
      const Conserved c = rotate_solve_rotate_back
	(rs_,
	 left,
	 reflect(left,p),
	 0,n,p);
      res.push_back
	(pair<size_t,Extensive>
	 (i,conserved_to_extensive
	  (c,left_cell)));
    }
  }
  return res;
}

namespace {

  Primitive boost(const Primitive& origin,
		  const Vector2D& v)
  {
    Primitive res = origin;
    res.Velocity += v;
    return res;
  }

  Conserved outflow_only
  (const RiemannSolver& rs,
   const Vector2D& rmp,
   const Edge& edge,
   const Primitive& cell,
   bool left_real)
  {
    const Vector2D& p = Parallel(edge);
    return left_real ?
      rotate_solve_rotate_back
      (rs,
       cell,
       cell.Velocity.y<0 ?
       reflect(cell,p) : cell,
       0,
       remove_parallel_component
       (edge.vertices.second-rmp,p),
       p) :
      rotate_solve_rotate_back
      (rs,
       cell.Velocity.y<0 ? 
       reflect(cell,p) : cell,
       cell,
       0,
       remove_parallel_component
       (rmp-edge.vertices.second,p),
       p);
  }
  
  Conserved support_riemann(const RiemannSolver& rs,
			    const Vector2D& rmp,
			    const Edge& edge,
			    const Primitive& cell,
			    const Vector2D& support,
			    bool left_real)
  {
    const Vector2D& p = Parallel(edge);
    Conserved res = left_real ?
      rotate_solve_rotate_back
      (rs,
       cell,
       boost(reflect(cell,p),support),
       0,
       remove_parallel_component
       (edge.vertices.second-rmp,p),
       p) :
      rotate_solve_rotate_back
      (rs,
       boost(reflect(cell,p),support),
       cell,
       0,
       remove_parallel_component
       (rmp-edge.vertices.second,p),
       p);
    res.Mass = 0;
    res.Energy = 0;
    return res;
  }

  Primitive gravinterpolate
  (const ComputationalCell& origin,
   const Vector2D& cm,
   const Vector2D& centroid,
   const Vector2D& acc,
   const EquationOfState& eos)
  {
    const double energy =
      eos.dp2e(origin.density,
	       origin.pressure,
	       origin.tracers);
    const double sound_speed =
      eos.dp2c(origin.density,
	       origin.pressure,
	       origin.tracers);
    return Primitive
      (origin.density,
       origin.pressure + origin.density*ScalarProd(acc,centroid-cm),
       origin.velocity,
       energy,
       sound_speed);
  }

  Conserved bulk_riemann
  (const RiemannSolver& rs,
   const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const Edge& edge,
   const CoreAtmosphereGravity::AccelerationCalculator& ac)
  {
    const Vector2D left_pos =
      tess.GetCellCM(edge.neighbors.first);
    const Vector2D right_pos =
      tess.GetCellCM(edge.neighbors.second);
    const Vector2D left_acc = ac(left_pos);
    const Vector2D right_acc = ac(right_pos);
    const Vector2D centroid =
      0.5*(edge.vertices.first+edge.vertices.second);
    const size_t left_index =
      static_cast<size_t>(edge.neighbors.first);
    const size_t right_index =
      static_cast<size_t>(edge.neighbors.second);
    const Primitive left =
      gravinterpolate
      (cells[left_index],
       left_pos,
       centroid,
       left_acc,
       eos);
    const Primitive right =
      gravinterpolate
      (cells[right_index],
       right_pos,
       centroid,
       right_acc,
       eos);
    const Vector2D p = Parallel(edge);
    const Vector2D n =
      tess.GetMeshPoint(edge.neighbors.second) -
      tess.GetMeshPoint(edge.neighbors.first);
    const double velocity = Projection
      (tess.CalcFaceVelocity
       (point_velocities.at(left_index),
	point_velocities.at(right_index),
	tess.GetCellCM(edge.neighbors.first),
	tess.GetCellCM(edge.neighbors.second),
	calc_centroid(edge)),n);
    return rotate_solve_rotate_back
      (rs,left,right,velocity,n,p);
  }

  double calc_radius_sqr(const Vector2D& p)
  {
    return ScalarProd(p,p);
  }

  bool point_below_edge
  (const Vector2D& p,
   const Edge& edge)
  {
    return calc_radius_sqr(p)<calc_radius_sqr(edge.vertices.first) &&
      calc_radius_sqr(p)<calc_radius_sqr(edge.vertices.second);
  }
}

const Conserved InnerBC::calcHydroFlux
(const Tessellation& tess,
 const vector<Vector2D>& point_velocities,
 const vector<ComputationalCell>& cells,
 const EquationOfState& eos,
 const size_t i,
 const CoreAtmosphereGravity::AccelerationCalculator& ac) const
{
  const Edge& edge = tess.GetEdge(static_cast<int>(i));
  const pair<bool,bool> flags
    (edge.neighbors.first>=0 && edge.neighbors.first<tess.GetPointNo(),
     edge.neighbors.second>=0 && edge.neighbors.second<tess.GetPointNo());
  assert(flags.first || flags.second);
  if(!flags.first)
    return support_riemann
      (rs_,
       tess.GetMeshPoint(edge.neighbors.second),
       edge,
       convert_to_primitive
       (cells.at
	(static_cast<size_t>(edge.neighbors.second)),
	eos),
       Vector2D(0,0),
       false);
  if(!flags.second)
    return support_riemann
      (rs_,
       tess.GetMeshPoint(edge.neighbors.first),
       edge,
       convert_to_primitive
       (cells.at
	(static_cast<size_t>(edge.neighbors.first)),
	eos),
       Vector2D(0,0),
       true);
  const size_t left_index =
    static_cast<size_t>(edge.neighbors.first);
  const size_t right_index =
    static_cast<size_t>(edge.neighbors.second);
  const ComputationalCell& left_cell = cells.at(left_index);
  const ComputationalCell& right_cell = cells.at(right_index);
  if(safe_retrieve(left_cell.stickers,ghost_)){
    if(safe_retrieve(right_cell.stickers,ghost_))
      return Conserved();
    return point_below_edge
      (tess.GetMeshPoint(edge.neighbors.second),edge) ?
      outflow_only
      (rs_,
       tess.GetMeshPoint(edge.neighbors.second),
       edge,
       convert_to_primitive(right_cell,eos),
       false) :
      support_riemann
      (rs_,
       tess.GetMeshPoint(edge.neighbors.second),
       edge,
       convert_to_primitive(right_cell,eos),
       Vector2D(0,0),
       false);
  }
  if(safe_retrieve(right_cell.stickers,ghost_)){
    return point_below_edge
      (tess.GetMeshPoint(edge.neighbors.first),edge) ?
      outflow_only
      (rs_,
       tess.GetMeshPoint(edge.neighbors.first),
       edge,
       convert_to_primitive(left_cell,eos),
       true) :
      support_riemann
      (rs_,
       tess.GetMeshPoint(edge.neighbors.first),
       edge,
       convert_to_primitive(left_cell,eos),
       Vector2D(0,0),
       true);
  }
  return bulk_riemann
    (rs_,
     tess,
     point_velocities,
     cells,
     eos,
     edge,
     ac);
}
