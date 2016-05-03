#include "temperature_cell_updater.hpp"

namespace
{
	void update_single(const FermiTable& eos, vector<Extensive>& extensives,
		const ComputationalCell& old,
		const CacheData& cd,
		const size_t index,
		ComputationalCell &res,
		size_t temperature_loc_,
		TracerStickerNames const& tracerstickernames)
	{
		Extensive& extensive = extensives[index];
		const double volume = cd.volumes[index];
		res.density = extensive.mass / volume;
		res.velocity = extensive.momentum / extensive.mass;
		const double energy = extensive.energy / extensive.mass -
			0.5*ScalarProd(res.velocity, res.velocity);
		res.stickers = old.stickers;
		for (size_t i = 0; i < extensive.tracers.size(); ++i)
			res.tracers[i] = extensive.tracers[i] / extensive.mass;
		try
		{
			double temperature = eos.de2t(res.density,energy,res.tracers,tracerstickernames.tracer_names);
			res.pressure = eos.dt2p(res.density,temperature,res.tracers,tracerstickernames.tracer_names);
			res.tracers[temperature_loc_] = temperature;
			extensive.tracers[temperature_loc_] = res.tracers[temperature_loc_] * extensive.mass;
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Cell index", static_cast<double>(index));
			eo.AddEntry("Cell mass", extensive.mass);
			eo.AddEntry("Cell x momentum", extensive.momentum.x);
			eo.AddEntry("Cell y momentum", extensive.momentum.y);
			eo.AddEntry("Cell energy", extensive.energy);
			throw;
		}
	}
}

TemperatureCellUpdater::TemperatureCellUpdater(string tracername,FermiTable const& fermi_eos)
: temperature_loc_(-1),tracername_(tracername),fermi_eos_(fermi_eos){}

vector<ComputationalCell> TemperatureCellUpdater::operator()
		(const Tessellation& tess,
			const PhysicalGeometry& pg,
			const EquationOfState& /*eos*/,
			vector<Extensive>& extensives,
			const vector<ComputationalCell>& old,
			const CacheData& cd,
			TracerStickerNames const& tracerstickernames) const
{
	if(temperature_loc_<0)
	{
		vector<string>::const_iterator it = binary_find(tracerstickernames.tracer_names.begin(),
			tracerstickernames.tracer_names.end(), tracername_);
		assert(it != tracerstickernames.tracer_names.end());
		temperature_loc_ = it - tracerstickernames.tracer_names.begin();
	}
	size_t N = static_cast<size_t>(tess.GetPointNo());
	vector<ComputationalCell> res(N, old[0]);
	string ghost_("ghost");
	for (size_t i = 0; i < N; ++i)
	{
		if(!safe_retrieve(old[i].stickers,tracerstickernames.sticker_names,ghost_))
			update_single(fermi_eos_, extensives, old[i], cd, i, res[i], temperature_loc_,tracerstickernames);
		else
			res[i]=old[i];
	}
	return res;
}