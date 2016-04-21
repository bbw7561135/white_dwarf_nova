#include <cmath>
#include <fstream>
#include "create_grid.hpp"
#include "rectangle_stretch.hpp"
#include "source/tessellation/right_rectangle.hpp"
#include "source/newtonian/test_2d/clip_grid.hpp"

using std::ofstream;
using std::endl;

vector<Vector2D> create_grid(const pair<Vector2D,Vector2D>& boundaries,
			     const double dq,
			     const double r_min)
{
  vector<Vector2D> res;
  //   const double r_min = 0.1*boundaries.first.y;
  const double r_max = boundaries.second.y;
  for(double r=r_min;
      r<r_max;
      r*=(1+dq)){
    for(double q=0;q<2*M_PI;q+=dq)
      res.push_back(Vector2D(r*cos(q),r*sin(q)));
  }
  //    res = rectangular_clip(res,rectangle_stretch(boundaries,0.99));
  res = clip_grid (RightRectangle(rectangle_stretch(boundaries,0.99)),res);
  /*
  ofstream f("mesh_points.txt");
  for(size_t i=0;i<res.size();++i)
    f << res.at(i).x << " " << res.at(i).y << endl;
  f.close();
  */
  return res;
}
