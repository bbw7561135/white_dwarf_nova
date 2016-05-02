

#ifndef TEMPERATURE_CELL_UPDATER_HPP
#define TEMPERATURE_CELL_UPDATER_HPP 1

#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "fermi_table.hpp"

using std::vector;

//! \brief Base class for cell update scheme
class TemperatureCellUpdater : public CellUpdater
{
public:
	TemperatureCellUpdater(string tracername,FermiTable const& fermi_eos);

	vector<ComputationalCell> operator()
		(const Tessellation& tess,
			const PhysicalGeometry& pg,
			const EquationOfState& eos,
			vector<Extensive>& extensives,
			const vector<ComputationalCell>& old,
			const CacheData& cd,
			TracerStickerNames const& tracerstickernames) const;
private:
	mutable int temperature_loc_;
	const string tracername_;
	FermiTable const& fermi_eos_;
};

#endif // TEMPERATURE_CELL_UPDATER_HPP
