#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
#include "plane.h"

#include <cassert>

//#include "container.h"
#include "geometry.h"
#include "planex.h"
#include "planey.h"
#include "planez.h"
// 
#pragma GCC diagnostic pop

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
  if (m_plane_z)
  {
    try
    {
      m_plane_z->GetFunctionA();
      m_plane_z->GetFunctionB();
      m_plane_z->GetFunctionC();
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
    max_error = std::max(max_error,m_plane_x->CalcMaxError(coordinat));
  }
  if (CanCalcY())
  {
    max_error = std::max(max_error,m_plane_y->CalcMaxError(coordinat));
  }
  if (CanCalcZ())
  {
    max_error = std::max(max_error,m_plane_z->CalcMaxError(coordinat));
  }
  return max_error;
}


ribi::Plane::Coordinats2D ribi::Plane::CalcProjection(
  const Coordinats3D& points
) const
{
  if (!m_plane_x && !m_plane_y && !m_plane_z)
  {
    throw std::logic_error("Plane::CalcProjection: cannot express any plane");
  }
  try { if (m_plane_x) { return m_plane_x->CalcProjection(points); }}
  catch (std::logic_error&) { /* OK, try next plane */ }

  try { if (m_plane_y) { return m_plane_y->CalcProjection(points); }}
  catch (std::logic_error&) { /* OK, try next plane */ }

  try { if (m_plane_z) { return m_plane_z->CalcProjection(points); }}
  catch (std::logic_error&) { /* OK, try next plane */ }

  // TRACE("ERROR");
  // TRACE("INITIAL POINTS");
  /* for (const auto& point: m_points)
  {
    std::stringstream s;
    s
      << "("
      << boost::geometry::get<0>(point) << ","
      << boost::geometry::get<1>(point) << ","
      << boost::geometry::get<2>(point)
      << ")"
    ;
    // TRACE(s.str());
  } */

  /* TRACE(points.size());
  {
    if (m_plane_x)
    {
      // try { TRACE(*m_plane_x); } catch(std::logic_error&) { TRACE("Failed m_plane_x"); }
      // try { m_plane_x->CalcProjection(points); } catch (std::logic_error&) { TRACE("Failed m_plane_x->CalcProjection"); }
    }
  }
  {
    if (m_plane_y)
    {
      // try { TRACE(*m_plane_y); } catch(std::logic_error&) { TRACE("Failed m_plane_y"); }
      // try { m_plane_y->CalcProjection(points); } catch (std::logic_error&) { TRACE("Failed m_plane_y->CalcProjection"); }
    }
  }
  {
    if (m_plane_z)
    {
      // try { TRACE(*m_plane_z); } catch(std::logic_error&) { TRACE("Failed m_plane_z"); }
      // try { m_plane_z->CalcProjection(points); } catch (std::logic_error&) { TRACE("Failed m_plane_z->CalcProjection"); }
    }
  } */
  /* for (const auto& point: points)
  {
    std::stringstream s;
    s
      << "("
      << boost::geometry::get<0>(point) << ","
      << boost::geometry::get<1>(point) << ","
      << boost::geometry::get<2>(point)
      << ")"
    ;
    // TRACE(s.str());
  } */
  assert(!"Should not get here");
  throw std::logic_error("Plane::CalcProjection: unexpected behavior");
}

ribi::Plane::Double ribi::Plane::CalcX(const Double& y, const Double& z) const
{
  if (!CanCalcX())
  {
    throw std::logic_error("Plane::CalcX: cannot express the plane as 'X = A*Y + B*Z + C'");
  }
  return m_plane_x->CalcX(y,z);
}

ribi::Plane::Double ribi::Plane::CalcY(const ribi::Plane::Double& x, const ribi::Plane::Double& z) const
{
  if (!CanCalcY())
  {
    throw std::logic_error("Plane::CalcY: cannot express the plane as 'Y = A*X + B*Y + C'");
  }
  return m_plane_y->CalcY(x,z);
}

ribi::Plane::Double ribi::Plane::CalcZ(const ribi::Plane::Double& x, const ribi::Plane::Double& y) const
{
  if (!CanCalcZ())
  {
    throw std::logic_error("Plane::CalcZ: cannot express the plane as 'Z = A*X + B*Y + C'");
  }
  return m_plane_z->CalcZ(x,y);
}

bool ribi::Plane::CanCalcX() const noexcept
{
  return m_plane_x.get();
}
bool ribi::Plane::CanCalcY() const noexcept
{
  return m_plane_y.get();
}

bool ribi::Plane::CanCalcZ() const noexcept
{
  return m_plane_z.get();
}

boost::shared_ptr<ribi::PlaneX> ribi::Plane::CreatePlaneX(
  const Coordinat3D& p1,
  const Coordinat3D& p2,
  const Coordinat3D& p3
) noexcept
{
  // const bool verbose{false};
  try
  {
    const boost::shared_ptr<PlaneX> p(
      boost::make_shared<PlaneX>(p1,p2,p3)
    );
    assert(p);
    return p;
  }
  catch (std::exception& e)
  {
    // if (verbose) { /* TRACE(e.what()); */ }
    return boost::shared_ptr<PlaneX>();
  }
}

boost::shared_ptr<ribi::PlaneY> ribi::Plane::CreatePlaneY(
  const Coordinat3D& p1,
  const Coordinat3D& p2,
  const Coordinat3D& p3
) noexcept
{
  try
  {
    const boost::shared_ptr<PlaneY> p
      = boost::make_shared<PlaneY>(p1,p2,p3);
    assert(p);
    return p;
  }
  catch (std::exception&)
  {
    return boost::shared_ptr<PlaneY>();
  }
}

boost::shared_ptr<ribi::PlaneZ> ribi::Plane::CreatePlaneZ(
  const Coordinat3D& p1,
  const Coordinat3D& p2,
  const Coordinat3D& p3
) noexcept
{
  try
  {
    const boost::shared_ptr<PlaneZ> p
      = boost::make_shared<PlaneZ>(p1,p2,p3);
    assert(p);
    return p;
  }
  catch (std::exception&)
  {
    return boost::shared_ptr<PlaneZ>();
  }
}

ribi::Plane::Doubles ribi::Plane::GetCoefficientsX() const
{
  if (!m_plane_x)
  {
    throw std::logic_error("Plane::GetCoefficientsX: cannot express the plane as 'X = A*Y + B*Z + C'");
  }
  return m_plane_x->GetCoefficients();
}

ribi::Plane::Doubles ribi::Plane::GetCoefficientsY() const
{
  if (!m_plane_y)
  {
    throw std::logic_error("Plane::GetCoefficientsY: cannot express the plane as 'Y = A*X + B*Z + C'");
  }
  return m_plane_y->GetCoefficients();
}

ribi::Plane::Doubles ribi::Plane::GetCoefficientsZ() const
{
  if (!m_plane_z)
  {
    throw std::logic_error("Plane::GetCoefficientsZ: cannot express the plane as 'Z = A*X + B*Y + C'");
  }
  return m_plane_z->GetCoefficients();
}


std::string ribi::Plane::GetVersion() noexcept
{
  return "1.9";
}

std::vector<std::string> ribi::Plane::GetVersionHistory() noexcept
{
  return {
    "2014-03-07: version 1.0: initial version",
    "2014-03-10: version 1.1: allow vertical planes",
    "2014-03-13: version 1.2: bug fixed",
    "2014-04-01: version 1.3: use of std::unique_ptr",
    "2014-06-13: version 1.4: added operator<<, ToStr calls operator<<, shortened time to compile",
    "2014-06-16: version 1.5: improved detection of planes that can be expressed in less than three dimensions"
    "2014-07-03: version 1.6: use of apfloat, improved accuracy",
    "2014-07-10: version 1.7: use of apfloat only",
    "2014-07-15: version 1.8: multiple bugfixes",
    "2014-08-02: version 1.9: use of stubs, to speed up testing"
  };
}

bool ribi::Plane::IsInPlane(const Coordinat3D& coordinat) const noexcept
{
  //return CalcError(coordinat) <= CalcMaxError(coordinat);;
  assert(m_plane_x || m_plane_y || m_plane_z);
  const bool is_in_plane {
       (m_plane_x && m_plane_x->IsInPlane(coordinat))
    || (m_plane_y && m_plane_y->IsInPlane(coordinat))
    || (m_plane_z && m_plane_z->IsInPlane(coordinat))
  };
  #ifndef NDEBUG
  const bool has_error_below_max = CalcError(coordinat) <= CalcMaxError(coordinat);
  /* if (is_in_plane != has_error_below_max)
  {
    // TRACE("ERROR");
    // TRACE(is_in_plane);
    // TRACE(has_error_below_max);
    // TRACE(CalcError(coordinat));
    // TRACE(CalcMaxError(coordinat));
    // TRACE("BREAK");
  } */
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
  if (plane.m_plane_x)
  {
    try { os << (*plane.m_plane_x); }
    catch (std::exception&) { os << "divnull"; }
  }
  else
  {
    os << "null";
  }
  os << ',';
  if (plane.m_plane_y)
  {
    try { os << (*plane.m_plane_y); }
    catch (std::exception&) { os << "divnull"; }
  }
  else
  {
    os << "null";
  }
  os << ',';
  if (plane.m_plane_z)
  {
    try { os << (*plane.m_plane_z); }
    catch (std::exception&) { os << "divnull"; }
  }
  else
  {
    os << "null";
  }
  return os;
}
