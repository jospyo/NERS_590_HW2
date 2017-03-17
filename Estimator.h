#ifndef _ESTIMATOR_HEADER_
#define _ESTIMATOR_HEADER_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <typeinfo>

#include "Particle.h"
#include "Material.h"
#include "Reaction.h"

class estimator {
  private:
    std::string estimator_name;
  protected:
    unsigned long long nhist;
  public:
     estimator( std::string label ) : estimator_name(label) {};
    ~estimator() {};

    virtual std::string name() final { return estimator_name; };

    virtual void score( particle* ) = 0;

    template< typename T >
    void score( particle*, T ) { assert(false); };

    void score( particle*, double, std::shared_ptr< material > ) {};

    virtual void endHistory()       = 0;
    virtual void report()           = 0;
};

class single_valued_estimator : public estimator {
  private:

  protected:
    double tally_hist, tally_sum, tally_squared;
  public:
     using estimator::score;

     single_valued_estimator(std::string label ) : estimator(label) { 
       nhist         = 0;
       tally_hist    = 0.0;   
       tally_sum     = 0.0; 
       tally_squared = 0.0;
     };
    ~single_valued_estimator() {};

     virtual void endHistory()    final { 
       nhist++;
       tally_sum     += tally_hist;
       tally_squared += tally_hist * tally_hist;
       tally_hist = 0.0; }

     virtual void score( particle* ) = 0;

     virtual void report() final {
       double mean = tally_sum / nhist;
       double var  = ( tally_squared / nhist - mean*mean ) / nhist;
       std::cout << name() << "   " << mean << "   " << std::sqrt( var ) / mean << std::endl;  
     };

};

class surface_current_estimator : public single_valued_estimator {
  private:

  public:
     surface_current_estimator( std::string label ) : single_valued_estimator(label) {};
    ~surface_current_estimator() {};

    void score( particle* );
};

class counting_estimator : public estimator {
  private:
    int count_hist;
    std::vector< double > tally;
  public:
     counting_estimator( std::string label ) : estimator(label) { count_hist = 0; };
    ~counting_estimator() {};

    void score( particle* );
    void endHistory();
    void report();
};

#endif
