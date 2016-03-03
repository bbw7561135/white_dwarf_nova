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
}

#endif // RICH_MPI

typedef pair<const ConditionActionSequence::Condition*,
	     const ConditionActionSequence::Action*> capp;

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
  proctess_(process_positions(outer_), outer_),
#endif // RICH_MPI
  tess_
  (ss ?
   ss->mesh_points :
   create_grid(outer_.getBoundary(),2e-3,0.9*id.radius_list.front()),
   outer_),
  eos_("eos_tab.coded",1,1,0,generate_atomic_properties()),
  rs_(),
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
   (capp(new IsBoundaryEdge, 
	 new RegularFlux(rs_)))
   ()),
  eu_(),
  cu_(),
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
   calc_init_cond(tess_,eos_,id,domain),
   eos_,
   point_motion_,
   evc_,
   force_,
   tsf_,
   fc_,
   eu_,
   cu_) 
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
