#ifndef _SOURCE_HEADER_
#define _SOURCE_HEADER_

#include <stack>
#include <memory>

#include "Point.h"
#include "Distribution.h"
#include "Particle.h"

class source {
  private:
    std::shared_ptr< distribution<point> > dist_pos;
    std::shared_ptr< distribution<point> > dist_dir;
  public:
     source( std::shared_ptr< distribution<point> > pos, std::shared_ptr< distribution<point> > dir )
       : dist_pos(pos), dist_dir(dir) {};
    ~source() {};

    std::stack<particle> sample();
};

#endif
