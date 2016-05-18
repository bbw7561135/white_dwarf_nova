#include "my_main_loop.hpp"
#include "temperature_appendix.hpp"
#include "volume_appendix.hpp"
#include "energy_appendix.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "source/misc/int2str.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include "write_cycle.hpp"
#include "source/newtonian/test_2d/multiple_diagnostics.hpp"
#include "nuclear_burn.hpp"
#include "atlas_support.hpp"
#include "filtered_conserved.hpp"
#include "multiple_manipulation.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#endif // RICH_MPI

using namespace simulation2d;


class ImprovedConstantTimeInterval: public Trigger
{
public:

  ImprovedConstantTimeInterval(double dt, double t_next=0): dt_(dt), t_next_(t_next),counter_(0) {}

  bool operator()(const hdsim& sim)
  {
	  double max_temp=0;
	  
	  if(sim.getTime()>t_next_)
	  {
		  t_next_ += dt_;
		  counter_++;
		  return true;
	  }
	  else
	  {
		  if(counter_<400)
		  {
		    TracerStickerNames tsn = sim.GetTracerStickerNames();
			size_t index = static_cast<size_t>(find(tsn.tracer_names.begin(), tsn.tracer_names.end(), string("Temperature"))-tsn.tracer_names.begin());
			vector<ComputationalCell> const& cells = sim.getAllCells();
			double max_temp = 0;
			size_t N = static_cast<size_t>(sim.getTessellation().GetPointNo());
			for (size_t i = 0; i < N; ++i)
				max_temp = max(max_temp, cells[i].tracers[index]);
#ifdef RICH_MPI
			double send=max_temp;
			MPI_Allreduce(&send,&max_temp,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif //RICH_MPI
			if (max_temp > 1e9)
				t_next_ = sim.getTime() + min(t_next_ - sim.getTime(),0.05*dt_);
			return false;
		  }
	  }
  }

private:
  const double dt_;
  mutable double t_next_;
  size_t counter_;
};

void my_main_loop(hdsim& sim, 
		  const FermiTable& eos, 
		  const BurnTime& tsf,
		  bool rerun)
{
 #ifdef RICH_MPI
	 int rank=0;
	 MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  write_snapshot_to_hdf5
    (sim,
     (rerun ? "rerun_initial" : "initial")+
     string("_")+int2str(rank)+".h5",
     vector<DiagnosticAppendix*>
     (1,new TemperatureAppendix(eos)));
#else
  write_snapshot_to_hdf5(sim,rerun ? "rerun_initial.h5" : "initial.h5",
			 vector<DiagnosticAppendix*>
			 (1,new TemperatureAppendix(eos)));
#endif //RICH_MPI
const double tf = 12;
//const double tf = 0.09;
  SafeTimeTermination term_cond(tf, 1e6);
  vector<DiagnosticFunction*> diag_list = VectorInitialiser<DiagnosticFunction*>()
    [new ConsecutiveSnapshots
     (new ImprovedConstantTimeInterval(tf/300),
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
  NuclearBurn* nb = new NuclearBurn
    (string("alpha_table"),
     string("ghost"),
     eos,
     string(rerun ? "rerun_burn_energy_history.txt" : "burn_energy_history.txt"));
  tsf.attach(nb);
  MultipleManipulation manip
    (VectorInitialiser<Manipulate*>
     (new AtlasSupport())
     (nb)
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
