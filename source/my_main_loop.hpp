#ifndef MY_MAIN_LOOP_HPP
#define MY_MAIN_LOOP_HPP 1

#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "fermi_table.hpp"
#include "burn_time.hpp"

void my_main_loop
(hdsim& sim,
 const FermiTable& eos,
 const BurnTime& tsf,
 bool rerun=false);

#endif // MY_MAIN_LOOP_HPP
