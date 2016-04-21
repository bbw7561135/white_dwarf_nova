#ifdef RICH_MPI
#include "source/mpi/MeshPointsMPI.hpp"
#include "source/misc/mesh_generator.hpp"
#endif
#include "sim_data.hpp"
#include "calc_bottom_area.hpp"

#ifdef RICH_MPI

namespace {
  vector<Vector2D> process_positions
  (const SquareBox& boundary)
  {
    int ws = 0;
    MPI_Comm_size(MPI_COMM_WORLD,&ws);
    const size_t ws2 = static_cast<size_t>(ws);
    const Vector2D lower_left = boundary.getBoundary().first;
    const Vector2D upper_right = boundary.getBoundary().second;
    return RandSquare
      (ws2,
       lower_left.x, 
       upper_right.x,
       lower_left.y, 
       upper_right.y);
  }

  size_t sort_point(Tessellation const& tess,
		    Vector2D const& point)
  {
    for(size_t i=0;i<static_cast<size_t>(tess.GetPointNo());++i){
      vector<Vector2D> vertices;
      ConvexHull(vertices,tess,static_cast<int>(i));
      if(PointInCell(vertices,point))
	return i;
    }
    assert(false && "Point does not belong in tessllation");
  }

  vector<Vector2D> distribute_grid
  (const Tessellation& proc_tess,
   const vector<Vector2D>& complete)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    vector<Vector2D> res;
    BOOST_FOREACH(const Vector2D& p, complete){
      if(static_cast<size_t>(rank)==
	 sort_point(proc_tess,p))
	res.push_back(p);
    }
    return res;
  }
}

#endif // RICH_MPI

typedef pair<const ConditionActionSequence::Condition*,
	     const ConditionActionSequence::Action*> capp;

typedef pair<const SimpleCellUpdater::Condition*,
	     const SimpleCellUpdater::Action*> scupp;

namespace{
  class BothSpecial: public ConditionActionSequence::Condition
  {
  public:

    explicit BothSpecial(const string& sticker_name):
      sticker_name_(sticker_name) {}

    virtual pair<bool, bool>
    operator()
    (const Edge& edge,
     const Tessellation& tess,
     const vector<ComputationalCell>& cells,
     const TracerStickerNames& tsn) const
    {
      const size_t neighbor_1_index = 
	static_cast<size_t>(edge.neighbors.first);
      const ComputationalCell cell_1 = cells.at(neighbor_1_index);
      const bool cond_1 = cell_1.stickers.front();

      const size_t neighbor_2_index =
	static_cast<size_t>(edge.neighbors.second);
      const ComputationalCell cell_2 = cells.at(neighbor_2_index);
      const bool cond_2 = cell_2.stickers.front();

      return pair<bool,bool>(cond_1 && cond_2, false);
    }

  private:
    const string sticker_name_;
  };

  class ZeroFlux: public ConditionActionSequence::Action
  {
  public:

    void operator()
    (const Edge& edge,
     const Tessellation& tess,
     const Vector2D& edge_velocity,
     const vector<ComputationalCell>& cells,
     const EquationOfState& eos,
     const bool aux,
     Extensive& res,
     double time,
     const TracerStickerNames& tsn) const
    {
      res.tracers = cells.at(0).tracers;
    }
  };
}

SimData::SimData(const InitialData& id,
		 const Units& u,
		 const CircularSection& domain,
		 const Snapshot* const ss):
  pg_(Vector2D(0,0), Vector2D(1,0)),
  outer_
  (Vector2D
   (1.2*domain.getRadii().second*cos(domain.getAngles().second),
    0.95*domain.getRadii().first),
   Vector2D
   (1.2*domain.getRadii().second*cos(domain.getAngles().first),
    1.05*domain.getRadii().second)),
#ifdef RICH_MPI
  proctess_
  (process_positions(outer_),
   outer_),
  tess_
  (proctess_,
   distribute_grid
   (proctess_,
    ss ?
    ss->mesh_points :
    create_grid(outer_.getBoundary(),2e-3,0.9*id.radius_list.front())),
   outer_),
#else
  tess_
  (ss ?
   ss->mesh_points :
   create_grid(outer_.getBoundary(),2e-3,0.9*id.radius_list.front()),
   outer_),
#endif // RICH_MPI
  eos_("eos_tab.coded",1,1,0,generate_atomic_properties()),
  rs_(),
  alt_point_motion_(),
  point_motion_(),
  evc_(),
  cag_
  (u.core_mass,
   linspace(id.radius_list.front(),id.radius_list.back(),100),
   u.gravitation_constant,
   domain.getAngles()),
  geom_force_(pg_.getAxis()),
  force_(VectorInitialiser<SourceTerm*>
	 (&cag_)
	 (&geom_force_)
	 ()),
  tsf_(0.3),
  fc_
  (VectorInitialiser<capp>
   (capp(new IsBoundaryEdge,
	 new RigidWallFlux(rs_)))
   (capp(new BothSpecial("ghost"),
	 new ZeroFlux))
   (capp(new RegularSpecialEdge("ghost"), 
	 new RigidWallFlux(rs_)))
   (capp(new IsBulkEdge,
	 new RegularFlux(rs_)))
   ()),
  eu_(),
  cu_
  (VectorInitialiser<scupp>
   (scupp(new HasSticker("ghost"), new SkipUpdate))()),
  init_cond_(calc_init_cond(tess_,eos_,id,domain)),
  #ifdef RICH_MPI
  procupdate_(outer_),
  #endif //RICH_MPI
  sim_
  (
#ifdef RICH_MPI
   proctess_,
#endif // RICH_MPI
   tess_,
   outer_,
   pg_,
   ss ?
   ss->cells :
   init_cond_.second,
   eos_,
   alt_point_motion_,
   evc_,
   force_,
   tsf_,
   fc_,
   eu_,
   cu_,
   init_cond_.first
   #ifdef RICH_MPI
   ,
   &procupdate_
#endif // RICH_MPI
   ) 
{
  if(ss){
    sim_.setStartTime(ss->time);
    delete ss;
  }
}

hdsim& SimData::getSim(void)
{
  return sim_;
}

const FermiTable& SimData::getEOS(void) const
{
  return eos_;
}
