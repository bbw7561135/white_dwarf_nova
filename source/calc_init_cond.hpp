#ifndef CALC_INIT_COND_HPP
#define CALC_INIT_COND_HPP 1

#include <vector>
#include "source/newtonian/two_dimensional/computational_cell_2d.hpp"
#include "source/tessellation/tessellation.hpp"
#include "fermi_table.hpp"
#include "initial_data.hpp"
#include "source/tessellation/shape_2d.hpp"

using std::vector;

pair<TracerStickerNames, vector<ComputationalCell> >
calc_init_cond
(const Tessellation& tess,
 const FermiTable& eos,
 const InitialData& id,
 const Shape2D& cd);

#endif // CALC_INIT_COND_HPP
