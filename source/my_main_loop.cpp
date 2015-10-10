#include "my_main_loop.hpp"
#include "temperature_appendix.hpp"
#include "volume_appendix.hpp"
#include "energy_appendix.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include "write_cycle.hpp"
#include "source/newtonian/test_2d/multiple_diagnostics.hpp"
#include "nuclear_burn.hpp"
#include "atlas_support.hpp"
#include "filtered_conserved.hpp"
#include "multiple_manipulation.hpp"

using namespace simulation2d;

void my_main_loop(hdsim& sim, const FermiTable& eos, bool rerun)
{
  write_snapshot_to_hdf5(sim,rerun ? "rerun_initial.h5" : "initial.h5",
			 vector<DiagnosticAppendix*>
			 (1,new TemperatureAppendix(eos)));
  const double tf = 10;
  SafeTimeTermination term_cond(tf, 1e6);
  vector<DiagnosticFunction*> diag_list = VectorInitialiser<DiagnosticFunction*>()
    [new ConsecutiveSnapshots
     (new ConstantTimeInterval(tf/1000),
      new Rubric(rerun ? "rerun_snapshot_" : "snapshot_",".h5"),
      VectorInitialiser<DiagnosticAppendix*>
      (new TemperatureAppendix(eos))
      (new EnergyAppendix(eos))
      (new VolumeAppendix())())]
    [new WriteTime(rerun ? "rerun_time.txt" : "time.txt")]
    [new WriteCycle(rerun ? "rerun_cycle.txt" : "cycle.txt")]
    [new FilteredConserved(rerun ? "rerun_total_conserved.txt" :"total_conserved.txt")]
    ();
  MultipleDiagnostics diag(diag_list);
  MultipleManipulation manip
    (VectorInitialiser<Manipulate*>
     (new AtlasSupport())
     (new NuclearBurn
      (string("alpha_table"),
       string("ghost"),
       eos,
       string(rerun ? "rerun_burn_energy_history.txt" : "burn_energy_history.txt")))
     ());
    main_loop(sim,
	    term_cond,
	    &hdsim::TimeAdvance2Heun,
	    &diag,
	    &manip);
    write_snapshot_to_hdf5(sim,rerun ? "rerun_final.h5" : "final.h5",
			   vector<DiagnosticAppendix*>
			   (1,new TemperatureAppendix(eos)));
}
