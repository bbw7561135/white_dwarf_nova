#include "rigid_sides.hpp"

namespace
{
	void ReverseNormalVelocity(ComputationalCell &cell, Edge const& edge, size_t index, Tessellation const& tess)
	{
		Vector2D normal;
		if (index == 1)
			normal = tess.GetMeshPoint(edge.neighbors.first) - tess.GetMeshPoint(edge.neighbors.second);
		else
			normal = tess.GetMeshPoint(edge.neighbors.second) - tess.GetMeshPoint(edge.neighbors.first);
		normal = normal / abs(normal);
		Vector2D toadd = -2 * ScalarProd(cell.velocity, normal)*normal;
		cell.velocity += toadd;
		return;
	}
	
	size_t IsBoundaryEdge(Edge const& edge, int npoints)
	{
		if (edge.neighbors.first >= npoints)
			return 1;
		if (edge.neighbors.second >= npoints)
			return 2;
		else
			return 0;
	}
	
	vector<std::pair<size_t, size_t> > GetOuterEdgesIndeces2(Tessellation const& tess,vector<ComputationalCell> const& cells)
	{
		vector<Edge> const& edges = tess.getAllEdges();
		vector<std::pair<size_t, size_t> > res;
		int npoints = tess.GetPointNo();
		for (size_t i = 0; i < edges.size(); ++i)
		{
			const size_t ghostindex = IsBoundaryEdge(edges[i], npoints);
			if(ghostindex==0) // deal with inside ghost points
			{
				if(cells[static_cast<size_t>(edges[i].neighbors.first)].stickers[0]&&
					!cells[static_cast<size_t>(edges[i].neighbors.second)].stickers[0])
				{
					if(tess.GetWidth(edges[i].neighbors.second)*0.2<edges[i].GetLength())
					{
						res.push_back(std::pair<size_t, size_t>(i, 1));
						continue;
					}
				}
				if(!cells[static_cast<size_t>(edges[i].neighbors.first)].stickers[0]&&
					cells[static_cast<size_t>(edges[i].neighbors.second)].stickers[0])
				{
					if(tess.GetWidth(edges[i].neighbors.first)*0.2<edges[i].GetLength())
					{
						res.push_back(std::pair<size_t, size_t>(i, 2));
						continue;
					}
				}
			}
#ifdef RICH_MPI
			if (ghostindex > 0)
			{
				bool not_inner_ghost = (ghostindex==1) ?  !cells.at(static_cast<size_t>(edges[i].neighbors.first)).stickers[0]
					: !cells.at(static_cast<size_t>(edges[i].neighbors.second)).stickers[0];
				if(!not_inner_ghost)
				{
					if(((ghostindex==1) ? tess.GetWidth(edges[i].neighbors.second) : tess.GetWidth(edges[i].neighbors.first))*0.2
						> edges[i].GetLength())
						continue;
				}
				if (tess.GetOriginalIndex(edges[i].neighbors.first) != tess.GetOriginalIndex(edges[i].neighbors.second) &&
					not_inner_ghost)
					continue;
			}
#endif
			if (ghostindex == 1)
				res.push_back(std::pair<size_t, size_t>(i, 1));
			if (ghostindex == 2)
				res.push_back(std::pair<size_t, size_t>(i, 2));
		}
		return res;
	}

}

boost::container::flat_map<size_t, ComputationalCell> RigidSides::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells, double /*time*/,TracerStickerNames const& /*tracerstickernames*/) const
{
	boost::container::flat_map<size_t, ComputationalCell> res;
	vector<std::pair<size_t, size_t> > ghosts = GetOuterEdgesIndeces2(tess,cells);
	for (size_t i = 0; i < ghosts.size(); ++i)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(ghosts[i].first));
		size_t ghost_index = ghosts[i].second == 1 ? static_cast<size_t> (edge.neighbors.first)
			: static_cast<size_t>(edge.neighbors.second);
		ComputationalCell ctemp = cells[static_cast<size_t>(ghosts[i].second == 2 ? edge.neighbors.first : edge.neighbors.second)];
		ReverseNormalVelocity(ctemp, edge, ghosts[i].second, tess);
		res.insert(std::pair<size_t, ComputationalCell>(ghost_index, ctemp));
	}
	return res;
}

Slope RigidSides::GetGhostGradient
(Tessellation const& tess,
	vector<ComputationalCell> const& cells,
	vector<Slope> const& /*gradients*/,
	size_t ghost_index, double /*time*/, Edge const& /*edge*/, TracerStickerNames const& /*tracerstickernames*/) const
{
	ComputationalCell cell(cells[static_cast<size_t>(tess.GetOriginalIndex(static_cast<int>(ghost_index)))]);
	cell.density = 0;
	cell.pressure = 0;
	cell.velocity = Vector2D(0, 0);
	cell.tracers = cells[0].tracers;
	size_t N = cell.tracers.size();
	for (size_t i = 0; i < N; ++i)
		cell.tracers[i] = 0;
	return Slope(cell, cell);
}
