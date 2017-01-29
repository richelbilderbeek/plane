#include "plane.h"


#include <boost/test/unit_test.hpp>

#include <boost/limits.hpp>
#include <boost/make_shared.hpp>

#include "container.h"
#include "geometry.h"
#include "planex.h"
#include "planey.h"
#include "planez.h"

using namespace ribi;
using Coordinat3D = ::ribi::Plane::Coordinat3D;

BOOST_AUTO_TEST_CASE(test_ribi_plane_construct_simple)
{
  using Plane = ribi::Plane;
  const boost::shared_ptr<Plane> plane{
    new Plane(
      Plane::Coordinat3D(0.0,0.0,0.0),
      Plane::Coordinat3D(0.0,1.0,0.0),
      Plane::Coordinat3D(1.0,0.0,0.0)
    )
  };
  const Plane plane2(
    Plane::Coordinat3D(0.0,0.0,0.0),
    Plane::Coordinat3D(0.0,1.0,0.0),
    Plane::Coordinat3D(1.0,0.0,0.0)
  );
}

BOOST_AUTO_TEST_CASE(test_ribi_plane_that_can_be_expressed_in_all_three_forms)
{
  using Coordinat3D = ::ribi::Plane::Coordinat3D;
  const Coordinat3D p1( 1.0, 2.0,3.0);
  const Coordinat3D p2( 4.0, 6.0,9.0);
  const Coordinat3D p3(12.0,11.0,9.0);
  const Plane p(p1,p2,p3);
  BOOST_CHECK(p.CanCalcX());
  BOOST_CHECK(p.CanCalcY());
  BOOST_CHECK(p.CanCalcZ());
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

BOOST_AUTO_TEST_CASE(test_ribi_plane_x_is_two)
{
  using Coordinat3D = ::ribi::Plane::Coordinat3D;
  const Coordinat3D p1(2.0, 2.0,3.0);
  const Coordinat3D p2(2.0, 6.0,9.0);
  const Coordinat3D p3(2.0,11.0,9.0);
  const Plane p(p1,p2,p3);
  //Can only calculate the X
  BOOST_CHECK( p.CanCalcX());
  BOOST_CHECK(!p.CanCalcY());
  BOOST_CHECK(!p.CanCalcZ());
  //Can project
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

BOOST_AUTO_TEST_CASE(test_ribi_plane_x_is_123)
{
  const Coordinat3D p1(123.0, 2.0,3.0);
  const Coordinat3D p2(123.0, 6.0,9.0);
  const Coordinat3D p3(123.0,11.0,9.0);
  const Plane p(p1,p2,p3);
  BOOST_CHECK(p.CanCalcX());
  BOOST_CHECK(!p.CanCalcY());
  BOOST_CHECK(!p.CanCalcZ());
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

BOOST_AUTO_TEST_CASE(test_ribi_plane_y_is_3)
{
  const Coordinat3D p1( 2.0, 3.0, 5.0);
  const Coordinat3D p2( 7.0, 3.0, 9.0);
  const Coordinat3D p3(11.0,3.0,13.0);
  const Plane p(p1,p2,p3);

  BOOST_CHECK(p.CanCalcY());
  BOOST_CHECK(!p.CanCalcX());
  BOOST_CHECK(!p.CanCalcZ());
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

BOOST_AUTO_TEST_CASE(test_ribi_plane_z_is_5)
{
  const Coordinat3D p1( 2.0, 3.0,5.0);
  const Coordinat3D p2( 7.0,11.0,5.0);
  const Coordinat3D p3(13.0,17.0,5.0);
  const Plane p(p1,p2,p3);
  BOOST_CHECK(p.CanCalcZ());
  BOOST_CHECK(!p.CanCalcX());
  BOOST_CHECK(!p.CanCalcY());
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

BOOST_AUTO_TEST_CASE(is_in_plane_for_z_is_0)
{
  //CanCalcZ and IsInPlane, Z = 0 plane, from 1.0 coordinat
  const Coordinat3D p1(0.0,0.0,0.0);
  const Coordinat3D p2(0.0,1.0,0.0);
  const Coordinat3D p3(1.0,0.0,0.0);
  const Plane p(p1,p2,p3);
  BOOST_CHECK(!p.CanCalcX());
  BOOST_CHECK(!p.CanCalcY());
  BOOST_CHECK( p.CanCalcZ());

  BOOST_CHECK( p.IsInPlane(Coordinat3D( 0.0, 0.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D(-1.0,-1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 1.0,-1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D(-1.0, 1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 1.0, 1.0,0.0)));
}

#ifdef CARE_ABOUT_SUCH_HIGH_PRECISION
BOOST_AUTO_TEST_CASE(is_in_plane_for_z_is_0_for_smallest_double)
{
  //CanCalcZ and IsInPlane, Z = 0 plane, from smallest possible coordinat");
  const double i = std::numeric_limits<double>::denorm_min();
  BOOST_CHECK(i > 0.0);
  const Coordinat3D p1(0.0,0.0,0.0);
  const Coordinat3D p2(0.0,  i,0.0);
  const Coordinat3D p3(  i,0.0,0.0);
  const Plane p(p1,p2,p3);
  BOOST_CHECK(!p.CanCalcX());
  BOOST_CHECK(!p.CanCalcY());
  BOOST_CHECK( p.CanCalcZ());

  BOOST_CHECK( p.IsInPlane(Coordinat3D( 0.0, 0.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D(-1.0,-1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 1.0,-1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D(-1.0, 1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 1.0, 1.0,0.0)));
}
#endif

#ifdef CARE_ABOUT_SUCH_HIGH_PRECISION
BOOST_AUTO_TEST_CASE(ribi_plane_0)
{
  const double i = std::numeric_limits<double>::max();
  BOOST_CHECK(i > 0.0);
  const Coordinat3D p1(0.0,0.0,0.0);
  const Coordinat3D p2(0.0,  i,0.0);
  const Coordinat3D p3(  i,0.0,0.0);
  const Plane p(p1,p2,p3);
  BOOST_CHECK(!p.CanCalcX());
  BOOST_CHECK(!p.CanCalcY());
  BOOST_CHECK( p.CanCalcZ());

  BOOST_CHECK( p.IsInPlane(Coordinat3D( 0.0, 0.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D(-1.0,-1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 1.0,-1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D(-1.0, 1.0,0.0)));
  BOOST_CHECK( p.IsInPlane(Coordinat3D( 1.0, 1.0,0.0)));
}
#endif // CARE_ABOUT_SUCH_HIGH_PRECISION

BOOST_AUTO_TEST_CASE(ribi_plane_1)
{
  //if (verbose) TRACE("CanCalcZ, Z = 0 plane, zooming in");
  //Sorted by difficulty
  const auto series = PlaneZ::GetTestSeries();
  for (const double i:series)
  {
    if (i == 0.0) continue;
    BOOST_CHECK(i != 0.0 && "Cannot express plane when all its coordinats are at origin");
    const Coordinat3D p1(0.0,0.0,0.0);
    const Coordinat3D p2(0.0,  i,0.0);
    const Coordinat3D p3(  i,0.0,0.0);
    const Plane p(p1,p2,p3);
    BOOST_CHECK(!p.CanCalcX());
    BOOST_CHECK(!p.CanCalcY());
    BOOST_CHECK( p.CanCalcZ());

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


BOOST_AUTO_TEST_CASE(ribi_plane_2)
{
  const bool show_warning{false};
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
  //if (verbose) TRACE("CanCalcZ, Z = z plane, varying height");
  const auto series = PlaneZ::GetTestSeries();
  for (const double z:series)
  {
    const double i = 1.0;
    const Coordinat3D p1(0.0,0.0,z);
    const Coordinat3D p2(0.0,  i,z);
    const Coordinat3D p3(  i,0.0,z);
    const Plane p(p1,p2,p3);
    BOOST_CHECK(!p.CanCalcX());
    BOOST_CHECK(!p.CanCalcY());
    BOOST_CHECK( p.CanCalcZ());

    for (const double j:series)
    {
      if (show_warning && !p.IsInPlane(Coordinat3D(j,j,z)))
      {
        std::stringstream s;
        s << "Warning: coordinat " << Geometry().ToStr(Coordinat3D(j,j,z))
          << " is determined not to be in a Plane that was created from points "
          << Geometry().ToStr(p1) << ", "
          << Geometry().ToStr(p2) << " and "
          << Geometry().ToStr(p3) << ". Error: "
          << p.CalcError(Coordinat3D(j,j,z))
          << ", max error: "
          << p.CalcMaxError(Coordinat3D(j,j,z))
          << " ("
          << (p.CalcError(Coordinat3D(j,j,z)) / p.CalcMaxError(Coordinat3D(j,j,z)))
          << "x)"
        ;
        std::cerr << s.str() << '\n';
      }
      if (!p.IsInPlane(Coordinat3D(j,j,z)))
      {
        const auto max_error = p.CalcMaxError(Coordinat3D(j,j,z));
        const auto error = p.CalcError(Coordinat3D(j,j,z));
        if (error == 0.0) continue;
        if (abs(1.0 - (max_error / error)) < 0.01)
        {
          //Allow another percent of freedom
          continue;
        }
      }
      BOOST_CHECK(p.IsInPlane(Coordinat3D(  j,  j,z)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D(  j, -j,z)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D( -j,  j,z)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D( -j, -j,z)));
    }
  }
}


BOOST_AUTO_TEST_CASE(ribi_plane_3)
{
  const auto series = PlaneZ::GetTestSeries();
  //if (verbose) TRACE("CanCalcZ, Z = z plane, zooming in");
  const double min = 10e+8;
  const double max = 10e-8;
  //for (const double z:series)
  for (double z=min; z<max; z*=10.0)
  {
    //for (const double i:series)
    for (double i=min; i<max; i*=10.0)
    {
      if (i == 0.0) continue;
      BOOST_CHECK(i != 0.0 && "Cannot express plane when all its coordinats are at origin");
      const Coordinat3D p1(0.0,0.0,z);
      const Coordinat3D p2(0.0,  i,z);
      const Coordinat3D p3(  i,0.0,z);
      const Plane p(p1,p2,p3);
      BOOST_CHECK(!p.CanCalcX());
      BOOST_CHECK(!p.CanCalcY());
      BOOST_CHECK( p.CanCalcZ());

      for (const double j:series)
      {
        if (!p.IsInPlane(Coordinat3D(j,j,z)))
        {
          std::stringstream s;
          s << "Warning: coordinat " << Geometry().ToStr(Coordinat3D(j,j,z))
            << " is determined not to be in a Plane that was created from points "
            << Geometry().ToStr(p1) << ", "
            << Geometry().ToStr(p2) << " and "
            << Geometry().ToStr(p3) << ". Error: "
            << p.CalcError(Coordinat3D(j,j,z))
            << ", max error: "
            << p.CalcMaxError(Coordinat3D(j,j,z))
            << " ("
            << (p.CalcError(Coordinat3D(j,j,z)) / p.CalcMaxError(Coordinat3D(j,j,z)))
            << "x)"
          ;
          std::cerr << s.str() << '\n';
          //continue;
        }
        BOOST_CHECK(p.IsInPlane(Coordinat3D(  j,  j,z)));
        BOOST_CHECK(p.IsInPlane(Coordinat3D(  j, -j,z)));
        BOOST_CHECK(p.IsInPlane(Coordinat3D( -j,  j,z)));
        BOOST_CHECK(p.IsInPlane(Coordinat3D( -j, -j,z)));
      }
    }
  }
}


BOOST_AUTO_TEST_CASE(ribi_plane_4)
{
  //IsInPlane for X=0 plane
  /*

     |####/
     B###C
     |##/
     |#/
     |/
  ---A------
    /|
   /#|
  /##|

  */

  //if (verbose) TRACE("IsInPlane, X = 0 plane, from 1.0 coordinats");
  const Coordinat3D p1(0.0,0.0,0.0);
  const Coordinat3D p2(0.0,0.0,1.0);
  const Coordinat3D p3(0.0,1.0,0.0);
  const Plane p(p1,p2,p3);
  BOOST_CHECK( p.CanCalcX());
  BOOST_CHECK(!p.CanCalcY());
  BOOST_CHECK(!p.CanCalcZ());

  for (double i = std::numeric_limits<double>::denorm_min();
    i < 1.0e-8; //std::numeric_limits<double>::max();
    i *= 10.0)
  {
    BOOST_CHECK(p.IsInPlane(Coordinat3D(0.0, i, i)));
    BOOST_CHECK(p.IsInPlane(Coordinat3D(0.0, i,-i)));
    BOOST_CHECK(p.IsInPlane(Coordinat3D(0.0,-i, i)));
    BOOST_CHECK(p.IsInPlane(Coordinat3D(0.0,-i,-i)));
  }
}

#ifdef CARE_ABOUT_SUCH_HIGH_PRECISION
BOOST_AUTO_TEST_CASE(ribi_plane_5)
{
  //if (verbose) TRACE("IsInPlane, X = 0 plane, from smallest possible coordinats");
  const double i = std::numeric_limits<double>::denorm_min();
  BOOST_CHECK(i > 0.0);
  const Coordinat3D p1(0.0,0.0,0.0);
  const Coordinat3D p2(0.0,0.0,i);
  const Coordinat3D p3(0.0,i,0.0);
  const Plane p(p1,p2,p3);
  BOOST_CHECK( p.CanCalcX());
  BOOST_CHECK(!p.CanCalcY());
  BOOST_CHECK(!p.CanCalcZ());
}
#endif // CARE_ABOUT_SUCH_HIGH_PRECISION


BOOST_AUTO_TEST_CASE(ribi_plane_6)
{
  //if (verbose) TRACE("IsInPlane, X = 0 plane, zooming in, #223");
  const auto series = PlaneZ::GetTestSeries();
  for (const double i:series)
  {
    if (i == 0.0) continue;
    BOOST_CHECK(i != 0.0 && "Cannot express plane when all its coordinats are at origin");
    const Coordinat3D p1(0.0,0.0,0.0);
    const Coordinat3D p2(0.0,0.0,  i);
    const Coordinat3D p3(0.0,  i,0.0);
    const Plane p(p1,p2,p3);

    for (const double j:series)
    {
      BOOST_CHECK(p.IsInPlane(Coordinat3D(0.0, j, j)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D(0.0, j,-j)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D(0.0,-j, j)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D(0.0,-j,-j)));
    }
  }
}

BOOST_AUTO_TEST_CASE(ribi_plane_7)
{
  //IsInPlane for Y=0 plane
  /*

     |####/
     B###/#
     |##/##
     |#/###
     |/####
  ---A---C-
    /|
   / |
  /  |

  */

  //if (verbose) TRACE("IsInPlane, Y = 0 plane, from 1.0 coordinats");
  const Coordinat3D p1(0.0,0.0,0.0);
  const Coordinat3D p2(0.0,0.0,1.0);
  const Coordinat3D p3(1.0,0.0,0.0);
  const Plane p(p1,p2,p3);
  BOOST_CHECK(!p.CanCalcX());
  BOOST_CHECK( p.CanCalcY());
  BOOST_CHECK(!p.CanCalcZ());

  for (double i = std::numeric_limits<double>::denorm_min();
    i < 1.0e8 /*std::numeric_limits<double>::max()*/;
    i *= 10.0
  )
  {
    BOOST_CHECK(p.IsInPlane(Coordinat3D( i,0.0, i)));
    BOOST_CHECK(p.IsInPlane(Coordinat3D( i,0.0,-i)));
    BOOST_CHECK(p.IsInPlane(Coordinat3D(-i,0.0, i)));
    BOOST_CHECK(p.IsInPlane(Coordinat3D(-i,0.0,-i)));
  }
}

#ifdef CARE_ABOUT_SUCH_HIGH_PRECISION
BOOST_AUTO_TEST_CASE(ribi_plane_8)
{
  //if (verbose) TRACE("IsInPlane, Y = 0 plane, from smallest possible coordinats");
  const double i = std::numeric_limits<double>::denorm_min();
  BOOST_CHECK(i > 0.0);
  const Coordinat3D p1(0.0,0.0,0.0);
  const Coordinat3D p2(0.0,0.0,  i);
  const Coordinat3D p3(  i,0.0,0.0);
  const Plane p(p1,p2,p3);
  BOOST_CHECK(!p.CanCalcX());
  BOOST_CHECK( p.CanCalcY());
  BOOST_CHECK(!p.CanCalcZ());
}
#endif // CARE_ABOUT_SUCH_HIGH_PRECISION

BOOST_AUTO_TEST_CASE(ribi_plane_9)
{
  //if (verbose) TRACE("IsInPlane, Y = 0 plane, zooming in, #223");
  const auto series = PlaneZ::GetTestSeries();
  for (const double i:series)
  {
    if (i == 0.0) continue;
    BOOST_CHECK(i != 0.0 && "Cannot express plane when all its coordinats are at origin");
    const Coordinat3D p1(0.0,0.0,0.0);
    const Coordinat3D p2(0.0,0.0,  i);
    const Coordinat3D p3(  i,0.0,0.0);
    const Plane p(p1,p2,p3);

    for (const double j:series)
    {
      BOOST_CHECK(p.IsInPlane(Coordinat3D( j,0.0, j)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D( j,0.0,-j)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D(-j,0.0, j)));
      BOOST_CHECK(p.IsInPlane(Coordinat3D(-j,0.0,-j)));
    }
  }
}


BOOST_AUTO_TEST_CASE(ribi_plane_10)
{
  //Create plane with different slopes
  /*
  |          /
  |         /
  |        /         __--D
  |       /      __--   /|
  |      /   __--      / |
  |     /__--         /  |
  |    B             /   +
  |   /         __--C   /
  |  /      __--    |  /
  | /   __--        | /
  |/__--            |/
  A-----------------+---------------

  A = p1 = Origin, fixed
  B = p2 = At y axis, fixed
  C = p3 = Above X axis, dependent on slope
  D = p4 = To be calculated
  */
  //if (verbose) TRACE("IsInPlane, Slope, Slope in Z direction");
  const auto series = PlaneZ::GetTestSeries();
  for (double slope = 1.0; slope > 1.0e-8; slope /= 10.0)
  {
    //if (verbose) { TRACE(slope); }
    const double slope_less = slope * 0.999999;
    const double slope_more = slope * 1.000001;

    const Coordinat3D p1(0.0,0.0,0.0);
    const Coordinat3D p2(0.0,1.0,0.0);
    const Coordinat3D p3(1.0,0.0,slope);
    const Plane p(p1,p2,p3);
    BOOST_CHECK( p.IsInPlane(Coordinat3D( 1.0, 1.0,slope)));

    if (slope_less < slope) //Not always true, when slope is very small
    {
      BOOST_CHECK(!p.IsInPlane(Coordinat3D( 1.0, 1.0,slope_less)));
    }
    if (slope_more > slope) //Not always true, when slope is very big
    {
      BOOST_CHECK(!p.IsInPlane(Coordinat3D( 1.0, 1.0,slope_more)));
    }
  }
}


BOOST_AUTO_TEST_CASE(ribi_plane_11)
{
  //Create plane with different slopes
  /*
  |          /
  |         /
  |        /__-D
  |     __--   |
  | __-- /     |
  B-    /      |
  |    /       |
  |   /        |
  |  /      __-C
  | /   __--  /
  |/__--     /
  A---------+-----------------------

  A = p1 = Origin, fixed
  B = p2 = At y axis, fixed
  C = p3 = Above X axis, dependent on slope
  D = p4 = To be calculated
  */
  //if (verbose) TRACE("IsInPlane, Slope, vertical plane");
  typedef std::pair<double,double> Pair;
  for (const Pair co:
    {
      std::make_pair(    1.0,    1.0),
      std::make_pair(   10.0,   10.0),
      std::make_pair(  100.0,  100.0),
      std::make_pair( 1000.0, 1000.0),
      std::make_pair(10000.0,10000.0),
      std::make_pair(1.0,1.000),
      std::make_pair(1.0,0.100),
      std::make_pair(1.0,0.010),
      std::make_pair(1.0,0.001),
    }
  )
  {
    const Coordinat3D p1(0.0,0.0,0.0);
    const Coordinat3D p2(0.0,0.0,1.0);
    const Coordinat3D p3(co.first,co.second,0.0);
    const Plane p(p1,p2,p3);
    BOOST_CHECK( p.IsInPlane(Coordinat3D(co.first,co.second,1.0)));
  }
}


BOOST_AUTO_TEST_CASE(ribi_plane_12)
{
  //if (verbose) TRACE("CalcProjection, from a crash in the program");
  const Coordinat3D p1( 1.0,-0.0,0.0);
  const Coordinat3D p2(-1.0, 0.0,0.0);
  const Coordinat3D p3( 1.0,-0.0,1.0);
  const Coordinat3D p4(-1.0, 0.0,1.0);
  const Plane p(p1,p2,p3);

  BOOST_CHECK(p.CanCalcY());
  BOOST_CHECK(!p.CanCalcX());
  BOOST_CHECK(!p.CanCalcZ());
  BOOST_CHECK(!p.CalcProjection( { p1,p2,p3,p4 } ).empty());
}


BOOST_AUTO_TEST_CASE(ribi_plane_13)
{
  //if (verbose) TRACE("IsInPlane, from #218");
  const double x1 = -0.5500000000000004884981308350688777863979339599609375;
  const double y1 = 2.000000000000000444089209850062616169452667236328125;
  const double z1 = 0; //left out the '.0' intentionally
  const double x2 = -3.78623595505618038004058689693920314311981201171875;
  const double y2 = 2; //left out the '.0' intentionally
  const double z2 = 0; //left out the '.0' intentionally
  const double x3 = -0.5500000000000004884981308350688777863979339599609375;
  const double y3 = 2.000000000000000444089209850062616169452667236328125;
  const double z3 = 10; //left out the '.0' intentionally
  const double x4 = -3.78623595505618038004058689693920314311981201171875;
  const double y4 = 2; //left out the '.0' intentionally
  const double z4 = 10; //left out the '.0' intentionally

  const Coordinat3D p1(x1,y1,z1);
  const Coordinat3D p2(x2,y2,z2);
  const Coordinat3D p3(x3,y3,z3);
  const Coordinat3D p4(x4,y4,z4);
  const Plane p(p1,p2,p3);
  try
  {
    const bool must_be_true = p.IsInPlane(p4);
    BOOST_CHECK(must_be_true);
  }
  catch (std::exception&)
  {
    std::stringstream s;
    s
      << std::setprecision(99)
      << p << '\n'
      << Container().ToStr(p.GetCoefficientsX()) << '\n'
      << Container().ToStr(p.GetCoefficientsY()) << '\n'
      << Container().ToStr(p.GetCoefficientsZ()) << '\n'
    ;
    std::cerr << s.str() << '\n';
  }
  BOOST_CHECK(p.IsInPlane(p4));
}


BOOST_AUTO_TEST_CASE(ribi_plane_14)
{
  //if (verbose) TRACE("IsInPlane, crashes with Plane v1.5");
  const double x1 = 0.0004035051226622692510832834944523028752882964909076690673828125;
  const double y1 = 0.00023296416881187433805568132161312178141088224947452545166015625;
  const double z1 = 0; //left out the '.0' intentionally

  const double x2 = 0.000403505141811931846741734464245610070065595209598541259765625;
  const double y2 = 0.00023296414405748076185791173298156309101614169776439666748046875;
  const double z2 = 0; //left out the '.0' intentionally

  const double x3 = 0.0004035051226622692510832834944523028752882964909076690673828125;
  const double y3 = 0.00023296416881187433805568132161312178141088224947452545166015625;
  const double z3 = 0.00025000000000000000520417042793042128323577344417572021484375;

  const double x4 = 0.000403505141811931846741734464245610070065595209598541259765625;
  const double y4 = 0.00023296414405748076185791173298156309101614169776439666748046875;
  const double z4 = 0.00025000000000000000520417042793042128323577344417572021484375;

  const Coordinat3D p1(x1,y1,z1);
  const Coordinat3D p2(x2,y2,z2);
  const Coordinat3D p3(x3,y3,z3);
  const Coordinat3D p4(x4,y4,z4);
  const Plane p(p1,p2,p3);

  BOOST_CHECK(p.IsInPlane(p4));
}

BOOST_AUTO_TEST_CASE(ribi_plane_15)
{
  //if (verbose) TRACE("IsInPlane, crashes with Plane v1.6");
  // TRACE '"ERROR"' line 392 in file '..\..\Classes\CppTriangleMesh\trianglemeshcellscreator.cpp': 'ERROR'
  // TRACE 's.str()' line 399 in file '..\..\Classes\CppTriangleMesh\trianglemeshcellscreator.cpp': '4
  // (-5,5,0) (index: 659) (-5,-0.9999999999999997779553950749686919152736663818359375,0) (index: 666) (-5,5,1) (index: 700) (-5,-0.9999999999999997779553950749686919152736663818359375,1) (index: 707) '
  // TRACE '"BREAK"' line 400 in file '..\..\Classes\CppTriangleMesh\trianglemeshcellscreator.cpp': 'BREAK'
  // Assertion failed!
  const double x1 = -5; //left out the '.0' intentionally
  const double y1 = 5; //left out the '.0' intentionally
  const double z1 = 0; //left out the '.0' intentionally

  const double x2 = -5; //left out the '.0' intentionally
  const double y2 = -0.9999999999999997779553950749686919152736663818359375;
  const double z2 = 0; //left out the '.0' intentionally

  const double x3 = -5; //left out the '.0' intentionally
  const double y3 = 5; //left out the '.0' intentionally
  const double z3 = 1; //left out the '.0' intentionally

  const double x4 = -5; //left out the '.0' intentionally
  const double y4 = -0.9999999999999997779553950749686919152736663818359375;
  const double z4 = 1; //left out the '.0' intentionally


  const Coordinat3D p1(x1,y1,z1);
  const Coordinat3D p2(x2,y2,z2);
  const Coordinat3D p3(x3,y3,z3);
  const Coordinat3D p4(x4,y4,z4);
  const Plane p(p1,p2,p3);
  try
  {
    BOOST_CHECK(p.IsInPlane(p4));
  }
  catch (std::exception&)
  {
    std::stringstream s;
    s
      << std::setprecision(99)
      << p << '\n'
      << Container().ToStr(p.GetCoefficientsX()) << '\n'
      << Container().ToStr(p.GetCoefficientsY()) << '\n'
      << Container().ToStr(p.GetCoefficientsZ()) << '\n'
    ;
    std::cerr << s.str() << '\n';
  }
  BOOST_CHECK(p.IsInPlane(p4));
}

