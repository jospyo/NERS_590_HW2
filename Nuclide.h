#ifndef _NUCLIDE_HEADER_
#define _NUCLIDE_HEADER_

#include <vector>
#include <string>
#include <memory>

#include "Reaction.h"

class nuclide {
  private:
    std::string nuclide_name;
    std::vector< std::shared_ptr< reaction > > rxn;
  public:
     nuclide( std::string label ) : nuclide_name(label) {};
    ~nuclide() {};

    std::string name() { return nuclide_name; }
    std::vector< std::shared_ptr< reaction > > getReactions() {return rxn;} ;

    void        addReaction( std::shared_ptr< reaction > );
    double      total_xs();

    std::shared_ptr< reaction > sample_reaction();
};


#endif
