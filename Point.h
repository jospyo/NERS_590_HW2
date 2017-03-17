#ifndef _POINT_HEADER_
#define _POINT_HEADER_

// point in 3d space
class point {
  public:
    double x, y, z;

     point( double a = 0.0, double b = 0.0, double c = 0.0 ) : x(a), y(b), z(c) {};
    ~point() {};

    void normalize();
};

// a ray in 3d space; has a point of origin and a direction unit vector
class ray {
  public:
    ray( point p, point d );
    ~ray() {};

    point pos;
    point dir;
};

#endif
