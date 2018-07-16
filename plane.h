#ifndef RIBI_PLANE_H
#define RIBI_PLANE_H

#include <vector>




#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#ifndef _WIN32
#include <boost/geometry/geometries/polygon.hpp>
#endif


#include "planex.h"
#include "planey.h"
#include "planez.h"

namespace ribi {

///Any 3D plane, even a single point
///Can be constructed from its equation and at least three 3D points
//A plane stores its coefficients in the following form:
// A.x + B.y + C.z = D
//Converting this to z being a function of x and y:
// -C.z =  A  .x + B  .y - D
//    z = -A/C.x - B/C.y + D/C
struct Plane
{
  typedef double Double;
  typedef boost::geometry::model::d2::point_xy<Double> Coordinat2D;
  typedef boost::geometry::model::point<Double,3,boost::geometry::cs::cartesian> Coordinat3D;
  typedef std::vector<Coordinat2D> Coordinats2D;
  typedef std::vector<Coordinat3D> Coordinats3D;
  typedef std::vector<Double> Doubles;

  ///Construct a Plane from three points
  /*

   |
   |   /
   +--+
   |\/|
   |/\|
 --+--+---------
  /|
 / |

  */
  explicit Plane(
    const Coordinat3D& p1,
    const Coordinat3D& p2,
    const Coordinat3D& p3
  ) noexcept;

  ///Get the 2D projection of these 3D points,
  ///Assumes these are in a Plane
  /*

  A: (0,0,1)                  A: (0,0)
  B: (1,0,0)                  B: (SQRT(2),0)
  C: (1,1,0)                  C: (SQRT(2),SQRT(2))

  |    /
  |   /                       |
  A-----C                     |  C
  |\/  /      -> becomes ->   | /|
  |/\ /                       |/ |
  +--B---                     A--B-----

  */
  Coordinats2D CalcProjection(const Coordinats3D& points) const;

  ///If the Plane can be expressed as X = A*Y + B*Z + C, return the X
  Double CalcX(const Double& y, const Double& z) const;

  ///If the Plane can be expressed as Y = A*X + B*Z + C, return the Y
  Double CalcY(const Double& y, const Double& z) const;

  ///If the Plane can be expressed as Z = A*X + B*Y + C, return the Z
  Double CalcZ(const Double& y, const Double& z) const;

  ///Can the Plane be expressed as X = A*Y + B*Z + C ?
  bool CanCalcX() const noexcept;

  ///Can the Plane be expressed as Y = A*X + B*Z + C ?
  bool CanCalcY() const noexcept;

  ///Can the Plane be expressed as Z = A*X + B*Y + C ?
  bool CanCalcZ() const noexcept;

  ///Calculates the error between plane and coordinat
  Double CalcError(const Coordinat3D& coordinat) const noexcept;

  ///Calculates the maximum allowed error for that coordinat for it to be in the plane
  Double CalcMaxError(const Coordinat3D& coordinat) const noexcept;

  ///If the Plane can be expressed as X = A*Y + B*Z + C, return the coefficients
  Doubles GetCoefficientsX() const;

  ///If the Plane can be expressed as Y = A*X + B*Z + C, return the coefficients
  Doubles GetCoefficientsY() const;

  ///If the Plane can be expressed as Z = A*X + B*Y + C, return the coefficients
  Doubles GetCoefficientsZ() const;

  ///Checks if the coordinat is in the plane
  bool IsInPlane(const Coordinat3D& coordinat) const noexcept;


  private:

  ///A non-horizontal plane; a plane that can be expressed as 'X(Y,Z) = A*Y + B*Z + C'
  const std::vector<PlaneX> m_plane_x;

  ///A non-horizontal plane; a plane that can be expressed as 'Y(X,Z) = A*X + B*Z + C'
  const std::vector<PlaneY> m_plane_y;

  ///A non-vertical plane; a plane that can be expressed as 'Z(X,Y) = A*X + B*Y + C'
  const std::vector<PlaneZ> m_plane_z;

  const Coordinats3D m_points;

  friend std::ostream& operator<<(std::ostream& os, const Plane& plane) noexcept;
};

///Creates a vector of size one with a PlaneX if it can
///be constructed from the coefficients. Creates
///an empty vector otherwise
std::vector<PlaneX> CreatePlaneX(const std::vector<double>& coefficients_x) noexcept;
std::vector<PlaneX> CreatePlaneX(
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p1,
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p2,
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p3
) noexcept;

///Creates a vector of size one with a PlaneY if it can
///be constructed from the coefficients. Creates
///an empty vector otherwise
std::vector<PlaneY> CreatePlaneY(const std::vector<double>& coefficients_y) noexcept;
std::vector<PlaneY> CreatePlaneY(
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p1,
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p2,
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p3
) noexcept;

///Creates a vector of size one with a PlaneZ if it can
///be constructed from the coefficients. Creates
///an empty vector otherwise
std::vector<PlaneZ> CreatePlaneZ(const std::vector<double>& coefficients_z) noexcept;
std::vector<PlaneZ> CreatePlaneZ(
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p1,
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p2,
  const boost::geometry::model::point<double,3,boost::geometry::cs::cartesian>& p3
) noexcept;

std::ostream& operator<<(std::ostream& os, const Plane& plane) noexcept;


} //~namespace ribi

#endif // RIBI_PLANE_H
