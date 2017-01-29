#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
#include "planez.h"

#include <cassert>

//#include "container.h"
#include "geometry.h"
// 
#pragma GCC diagnostic pop

ribi::PlaneZ::PlaneZ() noexcept
  : PlaneZ(
    Coordinat3D(0.0,0.0,0.0),
    Coordinat3D(1.0,0.0,0.0),
    Coordinat3D(0.0,1.0,0.0)
  )
{

}

ribi::PlaneZ::PlaneZ(
  const Doubles& coefficients
)
  : m_coefficients(coefficients)
{
  // const bool verbose{false};
  assert(GetCoefficients().size() == 4);

  if (m_coefficients[2] == 0.0)
  {
    // if (verbose) { /* TRACE(Container().ToStr(m_coefficients)); */}
    throw std::logic_error("PlaneZ (from coeffients) cannot be expressed in less than 3D space");
  }

  try
  {
    //TRACE(GetFunctionA());
    //TRACE(GetFunctionB());
    //TRACE(GetFunctionC());
    assert(GetFunctionA() == 0.0 || GetFunctionA() != 0.0);
    assert(GetFunctionB() == 0.0 || GetFunctionB() != 0.0);
    assert(GetFunctionC() == 0.0 || GetFunctionC() != 0.0);
  }
  catch (...)
  {
    assert(!"Should not get here");
  }
}

ribi::PlaneZ::PlaneZ(
  const Coordinat3D& p1,
  const Coordinat3D& p2,
  const Coordinat3D& p3
) : PlaneZ(CalcPlaneZ(p1,p2,p3))
{
  // const bool verbose{false};
  assert(GetCoefficients().size() == 4);

  if (m_coefficients[2] == 0.0)
  {
    // if (verbose) { /* TRACE(Container().ToStr(m_coefficients)); */}
    throw std::logic_error("Plane (from points) that can be expressed in less than 3D space");
  }
}

double ribi::PlaneZ::CalcError(const Coordinat3D& coordinat) const noexcept
{
  const double x = boost::geometry::get<0>(coordinat);
  const double y = boost::geometry::get<1>(coordinat);
  const double z = boost::geometry::get<2>(coordinat);
  const double expected{z};
  const double calculated{CalcZ(x,y)};
  const double error{std::abs(calculated - expected)};
  return error;
}

ribi::PlaneZ::Double ribi::PlaneZ::CalcMinErrorPerC() noexcept
{
  //min_error_per_c will be about 0.000000001
  //stub_value increases this jut a little, by a 0.000001%
  const double stub_value = 0.000000001 * 1.00000001;
  #define USE_STUB
  #ifdef USE_STUB
  return stub_value;
  #else //USE_STUB
  //PlaneZ calculates its own tolerance for errors, by measuring it
  static Double min_error_per_c = 0.0;
  if (min_error_per_c > 0.0) return min_error_per_c;

  //const double low = std::numeric_limits<double>::denorm_min();
  //const double high = std::numeric_limits<double>::max();
  const double low  = 1.0e-16;
  const double high = 1.0e+16;
  const double min_x = low;
  const double max_x = high;
  const double min_y = low;
  const double max_y = high;
  const double min_z = low;
  const double max_z = high;
  const double zero(0.0);

  for (double z = min_z; z < max_z; z*=10.0)
  {
    for (double y = min_y; y < max_y; y*=10.0)
    {
      for (double x = min_x; x < max_x; x*=10.0)
      {
        const Coordinat3D p1(0.0,0.0,z);
        const Coordinat3D p2(0.0,  y,z);
        const Coordinat3D p3(  x,0.0,z);
        const PlaneZ p(p1,p2,p3);
        for (const auto& p4: { p1,p2,p3 } )
        {
          const auto error = p.CalcError(p4);
          const auto error_per_c = error / p.GetFunctionC();
          assert(error_per_c >= zero);
          if (error_per_c > min_error_per_c)
          {
            min_error_per_c = error_per_c;
            //TRACE(min_error_per_c);
            //TRACE(x);
            //TRACE(y);
            //TRACE(z);
            //TRACE(p.GetFunctionC());
            //TRACE(double(min_error_per_c) / p.GetFunctionC());
            //std::stringstream s;
            //s << Geometry().ToStr(p4) << " " << min_error;
            //TRACE(s.str());
          }
        }
      }
    }
    //TRACE(min_error_per_c);
  }
  //TRACE(min_error_per_c);
  assert(min_error_per_c > zero);
  assert(min_error_per_c < stub_value);
  assert(min_error_per_c > 0.99 * stub_value);
  return min_error_per_c;
  #endif // USE_STUB
}

double ribi::PlaneZ::CalcMaxError(const Coordinat3D& /*coordinat*/) const noexcept
{
  assert(CalcMinErrorPerC() > double(0.0));
  const double max_error{std::abs(CalcMinErrorPerC() * GetFunctionC())};
  assert(max_error >= double(0.0));
  return max_error;
}

std::vector<double> ribi::PlaneZ::CalcPlaneZ(
  const Coordinat3D& p1,
  const Coordinat3D& p2,
  const Coordinat3D& p3
) noexcept
{
  
  const auto v(
    Geometry().CalcPlane(
      p1,
      p2,
      p3
    )
  );
  assert(v.size() == 4);
  return v;
}

ribi::PlaneZ::Coordinats2D ribi::PlaneZ::CalcProjection(
  const Coordinats3D& points
) const
{
  assert(points.size() >= 3);
  const double x_origin = 0.0;
  const double y_origin = 0.0;
  const double z_origin = CalcZ(x_origin,y_origin);

  Coordinats2D v;
  for (const auto& point: points)
  {
    const Double x(boost::geometry::get<0>(point));
    const Double y(boost::geometry::get<1>(point));
    const Double z(boost::geometry::get<2>(point));
    const Double dx =
      sqrt( //Apfloat does not add the std::
          ((x - x_origin) * (x - x_origin))
        + ((z - z_origin) * (z - z_origin))
      ) * (x - x_origin)
    ;
    const Double dy =
      sqrt( //Apfloat does not add the std::
          ((y - y_origin) * (y - y_origin))
        + ((z - z_origin) * (z - z_origin))
      ) * (y - y_origin)
    ;

    Coordinat2D point_xy(dx,dy);
    v.push_back(point_xy);
  }
  return v;
}

ribi::PlaneZ::Double ribi::PlaneZ::CalcZ(const Double& x, const Double& y) const
{
  // z = -A/C.x - B/C.y + D/C = (-A.x - B.y + D) / C
  // const bool verbose{false};
  const auto a = m_coefficients[0];
  const auto b = m_coefficients[1];
  const auto c = m_coefficients[2];
  const auto d = m_coefficients[3];
  if (c == 0.0)
  {
    throw std::logic_error("ribi::PlaneZ::CalcZ: cannot calculate Z of a vertical plane");
  }
  const double term1 = -a*x;
  const double term2 = -b*y;
  const double numerator = term1 + term2 + d;
  /* if (verbose)
  {
    // TRACE(numerator);
    // TRACE(c);
  } */
  const double result = numerator / c;
  return result;
}

double ribi::PlaneZ::GetFunctionA() const
{
  // const bool verbose{false};
  const auto coeff_a = m_coefficients[0];
  const auto coeff_c = m_coefficients[2];
  static const double zero(0.0);
  assert(coeff_c != zero);
  /* if (verbose)
  {
    // TRACE(coeff_c);
  } */
  const auto a = -coeff_a/coeff_c;
  return a;
}

double ribi::PlaneZ::GetFunctionB() const
{
  const auto coeff_b = m_coefficients[1];
  const auto coeff_c = m_coefficients[2];
  static const double zero(0.0);
  assert(coeff_c != zero);
  const auto b = -coeff_b/coeff_c;
  return b;
}

double ribi::PlaneZ::GetFunctionC() const
{
  const auto coeff_c = m_coefficients[2];
  const auto coeff_d = m_coefficients[3];
  static const double zero(0.0);
  assert(coeff_c != zero);

  try
  {
    const auto c = coeff_d/coeff_c;
    return c;
  }
  catch (...)
  {
    assert(!"Should not get here");
    throw;
  }
}

std::vector<double> ribi::PlaneZ::GetTestSeries() noexcept
{
  return
  {
     1.0,
    -1.0,
     std::numeric_limits<double>::epsilon(),
    -std::numeric_limits<double>::epsilon(),
     1.e8,
    -1.e8,
    0.0
    // std::numeric_limits<double>::denorm_min(),
    //-std::numeric_limits<double>::denorm_min(),
    // 1.e64,
    //-1.e64,
    //std::numeric_limits<double>::min(),
    //std::numeric_limits<double>::max()
  };
}

std::string ribi::PlaneZ::GetVersion() const noexcept
{
  return "1.6";
}

std::vector<std::string> ribi::PlaneZ::GetVersionHistory() const noexcept
{
  return {
    "2014-03-10: version 1.0: initial version, split off from Plane",
    "2014-03-10: version 1.1: bug fixed, only occurred at debugging",
    "2014-03-13: version 1.2: bug fixed",
    "2014-04-01: version 1.3: use of std::unique_ptr",
    "2014-07-03: version 1.4: use of apfloat",
    "2014-07-09: version 1.5: use double in interface only"
    "2014-07-10: version 1.6: use of apfloat only"
  };
}

bool ribi::PlaneZ::IsInPlane(const Coordinat3D& coordinat) const noexcept
{
  try
  {
    const double error = CalcError(coordinat);
    const double max_error = CalcMaxError(coordinat);
    return error <= max_error;
  }
  catch (std::exception& e)
  {
    // TRACE("ERROR");
    // TRACE(e.what());
    assert(!"Should not get here");
    throw;
  }
}

std::string ribi::PlaneZ::ToFunction() const
{
  std::stringstream s;
  s << (*this);
  return s.str();
}

std::ostream& ribi::operator<<(std::ostream& os, const PlaneZ& planez)
{
  try
  {
    os
      << "z=("
      << planez.GetFunctionA() << "*x) + ("
      << planez.GetFunctionB() << "*y) + "
      << planez.GetFunctionC()
    ;
  }
  catch (std::logic_error&)
  {
    const std::string error
      = "ribi::PlaneZ::ToFunction: cannot calculate Z of a vertical plane";
    throw std::logic_error(error.c_str());
  }
  return os;
}
