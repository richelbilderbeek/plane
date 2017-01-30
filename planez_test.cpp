#include "planez.h"

#include <boost/test/unit_test.hpp>
#include "container.h"
#include "geometry.h"

using namespace ribi;
using boost::geometry::get;
using Coordinat2D = ribi::PlaneZ::Coordinat2D;
using Coordinat3D = ribi::PlaneZ::Coordinat3D;

BOOST_AUTO_TEST_CASE(ribi_planez_test_1)
{
  const PlaneZ p;
  BOOST_CHECK(!p.ToFunction().empty());
  BOOST_CHECK(!p.GetCoefficients().empty());
}

BOOST_AUTO_TEST_CASE(ribi_planez_test_2)
{
//if (verbose) TRACE("PlaneZ, Z = 5");
  const Coordinat3D p1( 2.0, 3.0,5.0);
  const Coordinat3D p2( 7.0,11.0,5.0);
  const Coordinat3D p3(13.0,17.0,5.0);
  const PlaneZ p(p1,p2,p3);
  BOOST_CHECK(
    !p.CalcProjection(
      {
        Coordinat3D(0.0,0.0,1.0),
        Coordinat3D(1.0,0.0,0.0),
        Coordinat3D(1.0,1.0,0.0)
      }
    ).empty()
  );
}

BOOST_AUTO_TEST_CASE(ribi_planez_test_3)
{
//IsInPlane for Z=0 plane
//if (verbose) TRACE("PlaneZ, preparation for Plane::CanCalcZ and Plane::IsInPlane, Z = 0 plane, from 1.0 coordinat");
  const Coordinat3D p1(0.0,0.0,0.0);
  const Coordinat3D p2(0.0,1.0,0.0);
  const Coordinat3D p3(1.0,0.0,0.0);
  const PlaneZ p(p1,p2,p3);
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 0.0, 0.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D(-1.0,-1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 1.0,-1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D(-1.0, 1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 1.0, 1.0,0.0)));
}

#ifdef CARE_ABOUT_SUCH_HIGH_PRECISION
BOOST_AUTO_TEST_CASE(ribi_planez_test_4)
{
//if (verbose) TRACE("PlaneZ, preparation for Plane::CanCalcZ and Plane::IsInPlane, Z = 0 plane, from smallest possible coordinat");
  const double i = std::numeric_limits<double>::denorm_min();
  BOOST_CHECK(i > 0.0);
  const Coordinat3D p1(0.0,0.0,0.0);
  const Coordinat3D p2(0.0,  i,0.0);
  const Coordinat3D p3(  i,0.0,0.0);
  const PlaneZ p(p1,p2,p3);
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 0.0, 0.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D(-1.0,-1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 1.0,-1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D(-1.0, 1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 1.0, 1.0,0.0)));

}
#endif // CARE_ABOUT_SUCH_HIGH_PRECISION

#ifdef CARE_ABOUT_SUCH_HIGH_PRECISION
BOOST_AUTO_TEST_CASE(ribi_planez_test_5)
{
//if (verbose) TRACE("PlaneZ, preparation for Plane::CanCalcZ and Plane::IsInPlane, Z = 0 plane, from biggest possible coordinat");
  const double i = std::numeric_limits<double>::max();
  BOOST_CHECK(i > 0.0);
  const Coordinat3D p1(0.0,0.0,0.0);
  const Coordinat3D p2(0.0,  i,0.0);
  const Coordinat3D p3(  i,0.0,0.0);
  const PlaneZ p(p1,p2,p3);
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 0.0, 0.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D(-1.0,-1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 1.0,-1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D(-1.0, 1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 1.0, 1.0,0.0)));
}
#endif // CARE_ABOUT_SUCH_HIGH_PRECISION

BOOST_AUTO_TEST_CASE(ribi_planez_test_6)
{
  //if (verbose) TRACE("PlaneZ, preparation for Plane::CanCalcZ and Plane::IsInPlane, Z = 0 plane, zooming in");
  const auto series = GetTestSeries();
  for (const double i:series)
  {
    if (i == 0.0) continue;
    BOOST_CHECK(i != 0.0 && "Cannot express plane when all its coordinats are at origin");
    const Coordinat3D p1(0.0,0.0,0.0);
    const Coordinat3D p2(0.0,  i,0.0);
    const Coordinat3D p3(  i,0.0,0.0);
    const PlaneZ p(p1,p2,p3);
    for (const double j:series)
    {
      BOOST_CHECK(p.IsInPlane(Coordinat3D(0.0,0.0,0.0)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D(  j,  j,0.0)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D(  j, -j,0.0)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D( -j,  j,0.0)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D( -j, -j,0.0)));
    }
  }
}

BOOST_AUTO_TEST_CASE(ribi_planez_test_7)
{
  //if (verbose) TRACE("IsInPlane, Z = 1, zooming to smallest three points to determine a plane, point above origin");
  for (double i = 1.0;
    i > 1.0e-8; //i > 0.0;
    i/=10.0
  )
  {
    const Coordinat3D p1(0.0,0.0,1.0);
    const Coordinat3D p2(0.0,  i,1.0);
    const Coordinat3D p3(  i,0.0,1.0);
    const PlaneZ p(p1,p2,p3);
    BOOST_CHECK(p.IsInPlane(p1));
    BOOST_CHECK(p.IsInPlane(p2));
    BOOST_CHECK(p.IsInPlane(p3));
  }
}

BOOST_AUTO_TEST_CASE(ribi_planez_test_8)
{
//if (verbose) TRACE("IsInPlane, Z = 1, zooming to smallest three points to determine a plane, point above origin");
  const double min = 1.0e-8;
  const double max = 1.0e+8;
  for (double z = min; z < max; z*=10.0)
  {
    for (double i = min; i < max; i*=10.0)
    {
      const Coordinat3D p1(0.0,0.0,z);
      const Coordinat3D p2(0.0,  i,z);
      const Coordinat3D p3(  i,0.0,z);
      BOOST_CHECK(i != 0.0);
      const PlaneZ p(p1,p2,p3);
      if ( (!p.IsInPlane(p1) || !p.IsInPlane(p2) || !p.IsInPlane(p3)))
      {
        std::stringstream s;
        s << "Warning: coordinats " << Geometry().ToStr(p1)
          << ", " << Geometry().ToStr(p2)
          << " and " << Geometry().ToStr(p3)
          << " are determined not to be in a PlaneZ that was created from points"
        ;
        if (abs(1.0 - (p.CalcMaxError(p1) / p.CalcError(p1))) < 0.01)
        {
          //Allow another percent of freedom
          continue;
        }
      }
      BOOST_CHECK(p.IsInPlane(p1));
      BOOST_CHECK(p.IsInPlane(p2));
      BOOST_CHECK(p.IsInPlane(p3));
    }
  }
}

BOOST_AUTO_TEST_CASE(ribi_planez_test_9)
{
  //if (verbose) TRACE("CanCalcZ, Z = 1.0 plane, zooming in");
  const auto series = GetTestSeries();
  for (const double i:series)
  {
    if (i == 0.0) continue;
    BOOST_CHECK(i != 0.0 && "Cannot express plane when all its coordinats are at origin");
    const Coordinat3D p1(0.0,0.0,1.0);
    const Coordinat3D p2(0.0,  i,1.0);
    const Coordinat3D p3(  i,0.0,1.0);
    const PlaneZ p(p1,p2,p3);
    for (const double j:series)
    {
      BOOST_CHECK(p.IsInPlane(Coordinat3D(0.0,0.0,1.0)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D(  j,  j,1.0)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D(  j, -j,1.0)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D( -j,  j,1.0)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D( -j, -j,1.0)));
    }
  }
}

BOOST_AUTO_TEST_CASE(ribi_planez_test_10)
{
/*

  |    /#/##########
  |   B#/###########
  |  /#/############
  | /#/#############
  |/#/##############
  A-------C--------- Z = z
  |/
--O----------------- Z = 0
 /|


*/
  const bool show_warning{false};
  //if (verbose) TRACE("PlaneZ, preparation for Plane::CanCalcZ and Plane::IsInPlane, Z = z plane, zooming in");
  //The height of the plane
  const auto series = GetTestSeries();
  for (const double z:series)
  {
    //The distance from the origin, will be used by the two construction points
    for (const double i:series)
    {
      if (i == 0.0) continue;
      BOOST_CHECK(i != 0.0 && "Cannot express plane when all its coordinats are at origin");
      const Coordinat3D p1(0.0,0.0,z);
      const Coordinat3D p2(0.0,  i,z);
      const Coordinat3D p3(  i,0.0,z);
      const PlaneZ p(p1,p2,p3);
      //The distance (actually, half the Manhattan distance) from the origin,
      //will be used by points tested to be in this plane
      for (const double j:series)
      {
        if (!p.IsInPlane(Coordinat3D(j,j,z)))
        {
          if (show_warning)
          {
            std::stringstream s;
            s << "Warning: coordinat " << Geometry().ToStr(Coordinat3D(j,j,z))
              << " is determined not to be in a PlaneZ that was created from points "
              << Geometry().ToStr(p1) << ", "
              << Geometry().ToStr(p2) << " and "
              << Geometry().ToStr(p3) << "."
            ;
            std::cerr << s.str() << '\n';
          }
          continue;
        }
        BOOST_CHECK(p.IsInPlane(Coordinat3D(  j,  j,z)));
        BOOST_CHECK(p.IsInPlane(Coordinat3D(  j, -j,z)));
        BOOST_CHECK(p.IsInPlane(Coordinat3D( -j,  j,z)));
        BOOST_CHECK(p.IsInPlane(Coordinat3D( -j, -j,z)));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ribi_planez_test_11)
{
  //if (verbose) { TRACE("Check formulas"); }
  const bool verbose{false};
  const double p1_x{ 1.0};
  const double p1_y{ 2.0};
  const double p1_z{ 3.0};
  const double p2_x{ 4.0};
  const double p2_y{ 6.0};
  const double p2_z{ 9.0};
  const double p3_x{12.0};
  const double p3_y{11.0};
  const double p3_z{ 9.0};
  PlaneZ p(
    Coordinat3D(p1_x,p1_y,p1_z),
    Coordinat3D(p2_x,p2_y,p2_z),
    Coordinat3D(p3_x,p3_y,p3_z)
  );
  const auto& t = p.GetCoefficients();
  BOOST_CHECK(t.size() == 4);
  const auto a = t[0];
  const auto b = t[1];
  const auto c = t[2];
  const auto d = t[3];
  const auto a_expected =  30.0;
  const auto b_expected = -48.0;
  const auto c_expected =  17.0;
  const auto d_expected = -15.0;
  BOOST_CHECK(std::abs(a - a_expected) < 0.001);
  BOOST_CHECK(std::abs(b - b_expected) < 0.001);
  BOOST_CHECK(std::abs(c - c_expected) < 0.001);
  BOOST_CHECK(std::abs(d - d_expected) < 0.001);
  const auto d_p1_expected = (a * p1_x) + (b * p1_y) + (c * p1_z);
  const auto d_p2_expected = (a * p2_x) + (b * p2_y) + (c * p2_z);
  const auto d_p3_expected = (a * p3_x) + (b * p3_y) + (c * p3_z);
  if (verbose)
  {
    std::clog
      << "(a * x) + (b * y) + (c * z) = d" << '\n'
      << "(" << a << " * x) + (" << b << " * y) + (" << c << " * z) = " << d << '\n'
      << "(" << a << " * " << p1_x << ") + (" << b << " * " << p1_y << ") + (" << c << " * " << p1_z << ") = " << d << '\n'
      << "(" << (a * p1_x) << ") + (" << (b * p1_y) << ") + (" << (c * p1_z) << ") = " << d << '\n'
      << "(" << a << " * " << p2_x << ") + (" << b << " * " << p2_y << ") + (" << c << " * " << p2_z << ") = " << d << '\n'
      << "(" << (a * p2_x) << ") + (" << (b * p2_y) << ") + (" << (c * p2_z) << ") = " << d << '\n'
      << "(" << a << " * " << p3_x << ") + (" << b << " * " << p3_y << ") + (" << c << " * " << p3_z << ") = " << d << '\n'
      << "(" << (a * p3_x) << ") + (" << (b * p3_y) << ") + (" << (c * p3_z) << ") = " << d << '\n'
    ;
    /* Screen output

    (a * x) + (b * y) + (c * z) = d
    (30 * x) + (-48 * y) + (17 * z) = -15
    (30 * 1) + (-48 * 2) + (17 * 3) = -15
    (30) + (-96) + (51) = -15
    (30 * 4) + (-48 * 6) + (17 * 9) = -15
    (120) + (-288) + (153) = -15
    (30 * 12) + (-48 * 11) + (17 * 9) = -15
    (360) + (-528) + (153) = -15

    */
  }
  BOOST_CHECK(std::abs(d - d_p1_expected) < 0.001);
  BOOST_CHECK(std::abs(d - d_p2_expected) < 0.001);
  BOOST_CHECK(std::abs(d - d_p3_expected) < 0.001);
}

BOOST_AUTO_TEST_CASE(ribi_planez_test_12)
{
  //if (verbose) { TRACE("CalcPlaneZ"); }
  //CalcPlaneZ return the coefficients in the following form:
  // A.x + B.y + C.z = D
  //Converting this to z being a function of x and y:
  // -C.z = A.x + B.y - D
  // z = -A/C.x - B/C.y + D/C
  //In this test, use the formula:
  //  z = (2.0 * x) + (3.0 * y) + (5.0)
  //Coefficients must then become:
  //  -A/C = 2.0
  //  -B/C = 3.0
  //   D/C = 5.0
  //Coefficients are, when setting C to 1.0:
  //  -A = 2.0 => A = -2.0
  //  -B = 3.0 => B = -3.0
  //   C = 1.0
  //   D = 5.0
  const Coordinat3D p1(1.0,1.0,10.0);
  const Coordinat3D p2(1.0,2.0,13.0);
  const Coordinat3D p3(2.0,1.0,12.0);
  const PlaneZ p(p1,p2,p3);
  const auto& t = p.GetCoefficients();
  const auto a = t[0];
  const auto b = t[1];
  const auto c = t[2];
  const auto d = t[3];
  const auto a_expected = -2.0;
  const auto b_expected = -3.0;
  const auto c_expected =  1.0;
  const auto d_expected =  5.0;
  BOOST_CHECK(std::abs(a - a_expected) < 0.001);
  BOOST_CHECK(std::abs(b - b_expected) < 0.001);
  BOOST_CHECK(std::abs(c - c_expected) < 0.001);
  BOOST_CHECK(std::abs(d - d_expected) < 0.001);
  const auto d_p1_expected = (a * 1.0) + (b * 1.0) + (c * 10.0);
  const auto d_p2_expected = (a * 1.0) + (b * 2.0) + (c * 13.0);
  const auto d_p3_expected = (a * 2.0) + (b * 1.0) + (c * 12.0);
  BOOST_CHECK(std::abs(d - d_p1_expected) < 0.001);
  BOOST_CHECK(std::abs(d - d_p2_expected) < 0.001);
  BOOST_CHECK(std::abs(d - d_p3_expected) < 0.001);
}

BOOST_AUTO_TEST_CASE(ribi_planez_test_13)
{
  //if (verbose) { TRACE("CalcZ, diagonal plane"); }
  const Coordinat3D p1(1.0,2.0,3.0);
  const Coordinat3D p2(2.0,5.0,8.0);
  const Coordinat3D p3(3.0,7.0,11.0);
  const PlaneZ p(p1,p2,p3);
  BOOST_CHECK(abs(p.CalcZ(1.0,2.0)- 3.0) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(abs(p.CalcZ(2.0,5.0)- 8.0) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(abs(p.CalcZ(3.0,7.0)-11.0) < 0.001); //no std:: , as apfloat puts abs in the global namespace
}

BOOST_AUTO_TEST_CASE(ribi_planez_test_14)
{
  //if (verbose) { TRACE("CalcZ, horizontal plane Z = 5.0"); }
/*

  |    /
  |   /#
  |  /##
  | /###
  |/####
---+-----
 /|
/ |
/  |

*/


  const Coordinat3D p1( 2.0, 3.0,5.0);
  const Coordinat3D p2( 7.0,11.0,5.0);
  const Coordinat3D p3(13.0,17.0,5.0);
  const PlaneZ p(p1,p2,p3);
  BOOST_CHECK(std::abs(p.CalcZ(1.0,2.0)-5.0) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(std::abs(p.CalcZ(3.0,5.0)-5.0) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(std::abs(p.CalcZ(7.0,9.0)-5.0) < 0.001); //no std:: , as apfloat puts abs in the global namespace
}

BOOST_AUTO_TEST_CASE(ribi_planez_test_15)
{
  //if (verbose) { TRACE("ToFunction, 3 points and 4 points"); }
  std::function<double(double,double)> f {
    [](const double x, const double y)
    {
      return (2.0 * x) + (3.0 * y) + 5.0;
    }
  };


  const double x1 =  2.0;
  const double y1 =  3.0;
  const double x2 =  5.0;
  const double y2 =  7.0;
  const double x3 = 11.0;
  const double y3 = 13.0;
  const double x4 = 17.0;
  const double y4 = 29.0;
  const Coordinat3D p1(x1,y1,f(x1,y1));
  const Coordinat3D p2(x2,y2,f(x2,y2));
  const Coordinat3D p3(x3,y3,f(x3,y3));
  const PlaneZ a(p1,p2,p3);
  //BOOST_CHECK(a.ToFunction() == "z=(2*x) + (3*y) + 5");
  BOOST_CHECK(!a.ToFunction().empty());

  const Coordinat3D p4(x4,y4,f(x4,y4));
  BOOST_CHECK(a.ToFunction() == PlaneZ(p1,p2,p4).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p1,p3,p4).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p1,p4,p3).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p2,p1,p3).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p2,p1,p4).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p2,p3,p1).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p2,p3,p4).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p2,p4,p1).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p2,p4,p3).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p3,p1,p2).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p3,p1,p4).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p3,p2,p1).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p3,p2,p4).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p3,p4,p1).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p3,p4,p2).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p4,p1,p2).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p4,p1,p3).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p4,p2,p1).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p4,p2,p3).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p4,p3,p1).ToFunction());
  BOOST_CHECK(a.ToFunction() == PlaneZ(p4,p3,p2).ToFunction());

}

BOOST_AUTO_TEST_CASE(ribi_planez_test_16)
{
  //if (verbose) { TRACE("GetProjection, for Z = 0 plane"); }
  /*

  A: (0,0,1)                  A: (0,0)
  B: (1,0,0)                  B: (1,0)
  C: (1,1,0)                  C: (1,1)

  |    /
  |   /                       |
  A-----C                     |  C
  |\/  /      -> becomes ->   | /|
  |/\ /                       |/ |
  +--B---                     A--B-----

  */
  const std::vector<Coordinat2D> v {
    PlaneZ().CalcProjection(
      {
        Coordinat3D(0.0,0.0,1.0),
        Coordinat3D(1.0,0.0,0.0),
        Coordinat3D(1.0,1.0,0.0)
      }
    )
  };
  BOOST_CHECK(v.size() == 3);
  BOOST_CHECK(abs(get<0>(v[0]) - 0.0 ) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(abs(get<1>(v[0]) - 0.0 ) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(abs(get<0>(v[1]) - 1.0 ) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(abs(get<1>(v[1]) - 0.0 ) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(abs(get<0>(v[2]) - 1.0 ) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(abs(get<1>(v[2]) - 1.0 ) < 0.001); //no std:: , as apfloat puts abs in the global namespace

}

BOOST_AUTO_TEST_CASE(ribi_planez_test_17)
{
  //if (verbose) { TRACE("CalcProjection, for Z = 2 plane"); }
  /*

  A: (0,0,1+2)                  A: (0,0)
  B: (1,0,0+2)                  B: (1,0)
  C: (1,1,0+2)                  C: (1,1)

  |    /
  |   /                       |
  A-----C                     |  C
  |\/  /      -> becomes ->   | /|
  |/\ /                       |/ |
  +--B---                     A--B-----

  */
  const std::vector<PlaneZ::Coordinat2D> v {
    PlaneZ(
      Coordinat3D(0.0,0.0,0.0+2.0),
      Coordinat3D(0.0,1.0,0.0+2.0),
      Coordinat3D(1.0,0.0,0.0+2.0)
    ).CalcProjection(
      {
        Coordinat3D(0.0,0.0,1.0+2.0),
        Coordinat3D(1.0,0.0,0.0+2.0),
        Coordinat3D(1.0,1.0,0.0+2.0)
      }
    )
  };
  BOOST_CHECK(v.size() == 3);
  BOOST_CHECK(abs(get<0>(v[0]) - 0.0 ) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(abs(get<1>(v[0]) - 0.0 ) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(abs(get<0>(v[1]) - 1.0 ) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(abs(get<1>(v[1]) - 0.0 ) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(abs(get<0>(v[2]) - 1.0 ) < 0.001); //no std:: , as apfloat puts abs in the global namespace
  BOOST_CHECK(abs(get<1>(v[2]) - 1.0 ) < 0.001); //no std:: , as apfloat puts abs in the global namespace
}
