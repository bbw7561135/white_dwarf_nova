#ifndef SIM_DATA_HPP
#define SIM_DATA_HPP 1

#include "initial_data.hpp"
#include "units.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "fermi_table.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/source_terms/CenterGravity.hpp"
#include "monopole_self_gravity.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylindrical_complementary.hpp"
#include "source/newtonian/two_dimensional/source_terms/SeveralSources.hpp"
#include "lazy_extensive_updater.hpp"
#include "lazy_cell_updater.hpp"
#include "create_grid.hpp"
#include "generate_atomic_properties.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "circular_section.hpp"
#include "calc_init_cond.hpp"
#include "core_atmosphere_gravity.hpp"
#include "source/newtonian/two_dimensional/modular_flux_calculator.hpp"
//#include "source/newtonian/two_dimensional/ghost_point_generators/RigidWallGenerator.hpp"
#include "rigid_sides.hpp"
#include "safe_linear_gauss.hpp"
#include "reflective_ghost_throughout.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/condition_action_sequence_2.hpp"
#include "temperature_cell_updater.hpp"
#include "source/newtonian/two_dimensional/interpolations/LinearGaussImproved.hpp"
#include "conic_motion.hpp"
#include "burn_time.hpp"
#ifdef RICH_MPI
#include "source/mpi/SetLoad.hpp"
#endif // RICH_MPI

class SimData
{
public:

  SimData(const InitialData& id,
	  const Units& u,
	  const CircularSection& domain,
	  const Snapshot* const ss = 0);

  hdsim& getSim(void);

  const BurnTime& getTSF(void) const;

  const FermiTable& getEOS(void) const;

private:
  const CylindricalSymmetry pg_;
  const SquareBox outer_;
#ifdef RICH_MPI
VoronoiMesh proctess_;
#endif // RICH_MPI
  VoronoiMesh tess_;
  const FermiTable eos_;
  const Hllc rs_;
  Eulerian alt_point_motion_;
  ConicMotion point_motion_;
  const StationaryBox evc_;
  RigidSides ghosts_;
  LinearGaussImproved interp_;
  CoreAtmosphereGravity cag_;
  CylindricalComplementary geom_force_;
  SeveralSources force_;
  //  const SimpleCFL tsf_;
  const BurnTime tsf_;
  const ConditionActionSequence2 fc_;
  const SimpleExtensiveUpdater eu_;
  const TemperatureCellUpdater cu_;
  const pair<TracerStickerNames, vector<ComputationalCell> >
  init_cond_;
#ifdef RICH_MPI
  ConstNumberPerProc procupdate_;
#endif // RICH_MPI
  hdsim sim_;
};

#endif // SIM_DATA_HPP
