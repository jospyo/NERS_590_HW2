#ifndef _REACTION_HEADER_
#define _REACTION_HEADER_

#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include <stack>
#include <utility>

#include "Particle.h"
#include "Distribution.h"

class reaction {
  private:
    double rxn_xs;
  protected:
    std::string rxn_name;
  public:
     reaction( double x ) : rxn_xs(x) {};
    ~reaction() {};

    virtual std::string name() final { return rxn_name; };
    virtual double xs() final { return rxn_xs; };
    virtual void sample( particle* p, std::stack<particle>* bank ) = 0;
};

class capture_reaction : public reaction {
  private:
 
  public:
     capture_reaction( double x ) : reaction(x) { rxn_name = "capture"; };
    ~capture_reaction() {};

    void sample( particle* p, std::stack<particle>* bank );
};

class scatter_reaction : public reaction {
  private:
    std::shared_ptr< distribution<double> > scatter_dist; 
  public:
     scatter_reaction( double x, std::shared_ptr< distribution<double> > D ) :
       reaction(x), scatter_dist(D) { rxn_name = "scatter"; };
    ~scatter_reaction() {};

    void sample( particle* p, std::stack<particle>* bank );
};

class fission_reaction : public reaction {
  private:
    std::shared_ptr< distribution<int> >   multiplicity_dist; 
    std::shared_ptr< distribution<point> > isotropic;
  public:
     fission_reaction( double x, std::shared_ptr< distribution<int> > D ) :
       reaction(x), multiplicity_dist(D) { 
         rxn_name = "fission";
         isotropic = std::make_shared< isotropicDirection_distribution > ( "isotropic" ); 
       };
    ~fission_reaction() {};

    void sample( particle* p, std::stack<particle>* bank );
};

#endif
