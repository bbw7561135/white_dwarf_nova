#include <ctime>
#include "units.hpp"
#include "sim_data.hpp"
#include "my_main_loop.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include <fenv.h>
/*
#ifndef RICH_MPI
#include <mpi.h>
#endif // RICH_MPI
*/

using namespace std;

int main(int argc, char *argv[])
{
#ifdef RICH_MPI
  MPI_Init(NULL, NULL);
#endif // RICH_MPI
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  const clock_t begin = clock();
  const Units units;
  const InitialData id
    ("radius_list.txt",
     "density_list.txt",
     "temperature_list.txt",
     "velocity_list.txt");
  SimData sim_data(id,
		   units,
		   CircularSection(id.radius_mid.front(),
				   id.radius_mid.back(),
				   0.46*M_PI,
				   0.54*M_PI),
		   argc>1 ? new Snapshot(read_hdf5_snapshot(argv[1])) : 0);
  hdsim& sim = sim_data.getSim();
  my_main_loop(sim,sim_data.getEOS(),sim_data.getTSF(),argc>1);

  const clock_t end = clock();
  ofstream f("wall_time.txt");
  f << static_cast<double>(end-begin)/CLOCKS_PER_SEC << endl;
  f.close();

#ifdef RICH_MPI
  MPI_Finalize();
#endif // RICH_MPI

  return 0;
}
