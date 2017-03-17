#include <cmath>

#include "Point.h"

void point::normalize() {
  double norm = 1.0 / std::sqrt( x*x + y*y + z*z );
  x *= norm; y *= norm; z *= norm;
}

ray::ray( point p, point d ) : pos(p), dir(d) { dir.normalize(); }
