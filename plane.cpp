



#include "plane.h"

#include <cassert>

#include "geometry.h"
#include "planex.h"
#include "planey.h"
#include "planez.h"
// 


ribi::Plane::Plane(
  const Coordinat3D& p1,
  const Coordinat3D& p2,
  const Coordinat3D& p3
) noexcept
: m_plane_x(CreatePlaneX(p1,p2,p3)),
  m_plane_y(CreatePlaneY(p1,p2,p3)),
  m_plane_z(CreatePlaneZ(p1,p2,p3)),
  m_points( {p1,p2,p3} )
{
  if (CanCalcZ())
  {
    try
    {
      m_plane_z.front().GetFunctionA();
      m_plane_z.front().GetFunctionB();
      m_plane_z.front().GetFunctionC();
    }
    catch (...)
    {
      assert(!"Should not get here");
    }
  }

}

ribi::Plane::Double ribi::Plane::CalcError(const Coordinat3D& coordinat) const noexcept
{
  // const bool verbose{false};
  const double x{boost::geometry::get<0>(coordinat)};
  const double y{boost::geometry::get<1>(coordinat)};
  const double z{boost::geometry::get<2>(coordinat)};

  double min_error = std::numeric_limits<double>::max();

  //Absolute method
  if (CanCalcX())
  {
    const auto expected = x;
    const double calculated = CalcX(y,z);
    const double error{std::abs(calculated - expected)};
    min_error = std::min(error,min_error);
  }
  if (CanCalcY())
  {
    const auto expected = y;
    const auto calculated = CalcY(x,z);
    const double error{std::abs(calculated - expected)};
    min_error = std::min(error,min_error);
  }
  if (CanCalcZ())
  {
    const auto expected = z;
    const auto calculated = CalcZ(x,y);
    const double error{std::abs(calculated - expected)};
    min_error = std::min(error,min_error);
  }
  assert(min_error >= 0.0);
  return min_error;
}

ribi::Plane::Double ribi::Plane::CalcMaxError(const Coordinat3D& coordinat) const noexcept
{
  double max_error{std::numeric_limits<double>::denorm_min()};
  if (CanCalcX())
  {
    max_error = std::max(max_error,m_plane_x.front().CalcMaxError(coordinat));
  }
  if (CanCalcY())
  {
    max_error = std::max(max_error,m_plane_y.front().CalcMaxError(coordinat));
  }
  if (CanCalcZ())
  {
    max_error = std::max(max_error,m_plane_z.front().CalcMaxError(coordinat));
  }
  return max_error;
}


ribi::Plane::Coordinats2D ribi::Plane::CalcProjection(
  const Coordinats3D& points
) const
{
  if (!CanCalcX() && !CanCalcY() && !CanCalcZ())
  {
    throw std::logic_error("Plane::CalcProjection: cannot express any plane");
  }
  try { if (CanCalcX()) { return m_plane_x.front().CalcProjection(points); }}
  catch (std::logic_error&) { /* OK, try next plane */ }

  try { if (CanCalcY()) { return m_plane_y.front().CalcProjection(points); }}
  catch (std::logic_error&) { /* OK, try next plane */ }

  try { if (CanCalcZ()) { return m_plane_z.front().CalcProjection(points); }}
  catch (std::logic_error&) { /* OK, try next plane */ }


  assert(!"Should not get here");
  throw std::logic_error("Plane::CalcProjection: unexpected behavior");
}

ribi::Plane::Double ribi::Plane::CalcX(const Double& y, const Double& z) const
{
  if (!CanCalcX())
  {
    throw std::logic_error("Plane::CalcX: cannot express the plane as 'X = A*Y + B*Z + C'");
  }
  return m_plane_x.front().CalcX(y,z);
}

ribi::Plane::Double ribi::Plane::CalcY(const ribi::Plane::Double& x, const ribi::Plane::Double& z) const
{
  if (!CanCalcY())
  {
    throw std::logic_error("Plane::CalcY: cannot express the plane as 'Y = A*X + B*Y + C'");
  }
  return m_plane_y.front().CalcY(x,z);
}

ribi::Plane::Double ribi::Plane::CalcZ(const ribi::Plane::Double& x, const ribi::Plane::Double& y) const
{
  if (!CanCalcZ())
  {
    throw std::logic_error("Plane::CalcZ: cannot express the plane as 'Z = A*X + B*Y + C'");
  }
  return m_plane_z.front().CalcZ(x,y);
}

bool ribi::Plane::CanCalcX() const noexcept
{
  return !m_plane_x.empty();
}
bool ribi::Plane::CanCalcY() const noexcept
{
  return !m_plane_y.empty();
}

bool ribi::Plane::CanCalcZ() const noexcept
{
  return !m_plane_z.empty();
}

std::vector<ribi::PlaneX> ribi::CreatePlaneX(
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p1,
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p2,
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p3
) noexcept
{
  std::vector<PlaneX> p;
  try
  {
    p.push_back(PlaneX(p1,p2,p3));
  }
  catch (std::exception&) {} //!OCLINT this is ok
  return p;
}

std::vector<ribi::PlaneY> ribi::CreatePlaneY(
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p1,
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p2,
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p3
) noexcept
{
  std::vector<PlaneY> p;
  try
  {
    p.push_back(PlaneY(p1,p2,p3));
  }
  catch (std::exception&) {} //!OCLINT this is ok
  return p;
}

std::vector<ribi::PlaneZ> ribi::CreatePlaneZ(
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p1,
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p2,
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p3
) noexcept
{
  std::vector<PlaneZ> p;
  try
  {
    p.push_back(PlaneZ(p1,p2,p3));
  }
  catch (std::exception&) {} //!OCLINT this is ok
  return p;
}

ribi::Plane::Doubles ribi::Plane::GetCoefficientsX() const
{
  if (!CanCalcX())
  {
    throw std::logic_error("Plane::GetCoefficientsX: cannot express the plane as 'X = A*Y + B*Z + C'");
  }
  return m_plane_x.front().GetCoefficients();
}

ribi::Plane::Doubles ribi::Plane::GetCoefficientsY() const
{
  if (!CanCalcY())
  {
    throw std::logic_error("Plane::GetCoefficientsY: cannot express the plane as 'Y = A*X + B*Z + C'");
  }
  return m_plane_y.front().GetCoefficients();
}

ribi::Plane::Doubles ribi::Plane::GetCoefficientsZ() const
{
  if (!CanCalcZ())
  {
    throw std::logic_error("Plane::GetCoefficientsZ: cannot express the plane as 'Z = A*X + B*Y + C'");
  }
  return m_plane_z.front().GetCoefficients();
}

bool ribi::Plane::IsInPlane(const Coordinat3D& coordinat) const noexcept
{
  //return CalcError(coordinat) <= CalcMaxError(coordinat);
  assert(CanCalcX() || CanCalcY() || CanCalcZ());
  const bool is_in_plane {
       (CanCalcX() && m_plane_x.front().IsInPlane(coordinat))
    || (CanCalcY() && m_plane_y.front().IsInPlane(coordinat))
    || (CanCalcZ() && m_plane_z.front().IsInPlane(coordinat))
  };
  #ifndef NDEBUG
  const bool has_error_below_max = CalcError(coordinat) <= CalcMaxError(coordinat);
  assert(is_in_plane == has_error_below_max);
  #endif
  return is_in_plane;
}

std::ostream& ribi::operator<<(std::ostream& os, const Plane& plane) noexcept
{
  os << '(';
  const auto n_points = static_cast<int>(plane.m_points.size());
  for (/* const */ auto i=0; i!=n_points; ++i)
  {
    os << Geometry().ToStr(plane.m_points[i]);
    os << (i != n_points - 1 ? ',' : ')');
  }
  os << ',';
  if (!plane.m_plane_x.empty())
  {
    try { os << plane.m_plane_x.front(); }
    catch (std::exception&) { os << "divnull"; }
  }
  else
  {
    os << "null";
  }
  os << ',';
  if (!plane.m_plane_y.empty())
  {
    try { os << plane.m_plane_y.front(); }
    catch (std::exception&) { os << "divnull"; }
  }
  else
  {
    os << "null";
  }
  os << ',';
  if (!plane.m_plane_z.empty())
  {
    try { os << plane.m_plane_z.front(); }
    catch (std::exception&) { os << "divnull"; }
  }
  else
  {
    os << "null";
  }
  return os;
}
