#include "sim_data.hpp"
#include "calc_bottom_area.hpp"

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
  tess_
  (ss ?
   ss->mesh_points :
   create_grid(outer_.getBoundary(),2e-3,0.9*id.radius_list.front()),
   outer_),
  eos_("eos_tab.coded",1,1,0,generate_atomic_properties()),
  rs_(),
  point_motion_(),
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
  gpg_("ghost"),
  sr_
  (eos_,
   gpg_,
   true,
   0.2,
   0.5,
   0.7,
   VectorInitialiser<string>
   ("He4")
   ("C12")
   ("O16")
   ("Ne20")
   ("Mg24")
   ("Si28")
   ("S32")
   ("Ar36")
   ("Ca40")
   ("Ti44")
   ("Cr48")
   ("Fe52")
   ("Ni56")()),
  hbc_(rs_,"ghost",eos_),
  fc_(gpg_,sr_,rs_,hbc_),
  eu_(),
  cu_(),
  sim_(tess_,
       outer_,
       pg_,
       ss ?
       ss->cells :
       calc_init_cond(tess_,eos_,id,domain),
       eos_,
       point_motion_,
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
