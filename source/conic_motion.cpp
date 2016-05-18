#include "conic_motion.hpp"

ConicMotion::ConicMotion(const CircularSection& domain,const EquationOfState& eos):circ_section_(domain),bpm_(Lagrangian()),
pm_(RoundCells(bpm_,eos)){}

namespace
{
	void SlowStickers(const Tessellation& tess,const vector<ComputationalCell>& cells,const CircularSection& circ_section_,
		vector<Vector2D> &res)
	{
		size_t Nloop = res.size();
		vector<int> neigh;
		Vector2D norm_right = Vector2D(-sin(circ_section_.getAngles().second),cos(circ_section_.getAngles().second));
		Vector2D norm_left = Vector2D(sin(circ_section_.getAngles().first),-cos(circ_section_.getAngles().first));
		for(size_t i=0;i<Nloop;++i)
		{
			if(cells[i].stickers[0])
			{
				res[i]=Vector2D(0,0);
				continue;
			}
			tess.GetNeighbors(static_cast<int>(i),neigh);
			for(size_t j=0;j<neigh.size();++j)
			{
				if(cells[static_cast<size_t>(neigh[j])].stickers[0])
				{
					res[i]=Vector2D(0,0);
					continue;
				}
			}
			double R = tess.GetWidth(static_cast<int>(i));
			double r = abs(tess.GetMeshPoint(static_cast<int>(i)));
			double theta = atan2(tess.GetMeshPoint(static_cast<int>(i)).y,tess.GetMeshPoint(static_cast<int>(i)).x);
			if(r*(circ_section_.getAngles().second - theta) < 2 * R && ScalarProd(res[i],norm_right) < 0)
			{
				res[i] = res[i] * 0.2;
				continue;
			}
			if(r*(theta - circ_section_.getAngles().first) < 2 * R && ScalarProd(res[i],norm_left) < 0)
			{
				res[i] = res[i] * 0.2;
				continue;
			}
		}
	}
	
}


vector<Vector2D> ConicMotion::operator()(const Tessellation& tess,const vector<ComputationalCell>& cells,double time,
TracerStickerNames const& tracerstickernames) const
{
	vector<Vector2D> res = pm_(tess,cells,time,tracerstickernames);
	return res;
}

vector<Vector2D> ConicMotion::ApplyFix(Tessellation const& tess, vector<ComputationalCell> const& cells, double time,
double dt, vector<Vector2D> const& velocities, TracerStickerNames const& tracerstickernames)const
{
	vector<Vector2D> res = pm_.ApplyFix(tess,cells,time,dt,velocities,tracerstickernames);
	SlowStickers(tess,cells,circ_section_,res);
	return res;
}

