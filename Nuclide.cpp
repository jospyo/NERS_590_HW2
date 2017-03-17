#include <vector>
#include <memory>
#include <cassert>

#include "Random.h"
#include "Nuclide.h"

// add a new reaction to the current nuclide
void nuclide::addReaction( std::shared_ptr< reaction > R ) { rxn.push_back( R ); }

// return the total microscopic cross section
double nuclide::total_xs() {
  double xs = 0.0;
  for ( auto r : rxn ) { xs += r->xs(); }
  return xs;
}

// randomly sample a reaction type from this nuclide
std::shared_ptr< reaction > nuclide::sample_reaction() {
  double u = total_xs() * Urand();
  double s = 0.0;
  for ( auto r : rxn ) {
    s += r->xs();
    if ( s > u ) { return r; }
  }
  assert( false ); // should never reach here
  return nullptr;
}
