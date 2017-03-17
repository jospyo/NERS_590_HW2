#include <vector>
#include <utility>
#include <memory>
#include <cassert>

#include "Random.h"
#include "Particle.h"
#include "Reaction.h"
#include "Nuclide.h"
#include "Material.h"

// add a nuclide and atom fraction to the material
void material::addNuclide( std::shared_ptr< nuclide > N, double frac ) { 
  nuclides.push_back( std::make_pair( N, frac ) ); 
}

// private utility function that returns sum of atomic fraction * microscopic total xs
// multiply this by atomic density to get macroscopic cross section
// (this function is useful to accelerate sampling the nuclide)
double material::micro_xs() {
  double xs = 0.0;
  for ( auto n : nuclides ) { 
    // first is pointer to nuclide, second is atomic fraction
    xs += n.first->total_xs() * n.second;
  }
  return xs;
}

// return the macroscopic cross section
double material::macro_xs() {
  return atom_density() * micro_xs();
}

// randomly sample a nuclide based on total cross sections and atomic fractions
std::shared_ptr< nuclide > material::sample_nuclide() {
  double u = micro_xs() * Urand();
  double s = 0.0;
  for ( auto n : nuclides ) {
    // first is pointer to nuclide, second is atomic fraction
    s += n.first->total_xs() * n.second;
    if ( s > u ) { return n.first; }
  }
  assert( false ); // should never reach here
  return nullptr;
}

// function that samples an entire collision: sample nuclide, then its reaction, 
// and finally process that reaction with input pointers to the working particle p
// and the particle bank
void material::sample_collision( particle* p, std::stack<particle>* bank ) {
  // first sample nuclide
  std::shared_ptr< nuclide >  N = sample_nuclide();

  // now get the reaction
  std::shared_ptr< reaction > R = N->sample_reaction();

  // finally process the reaction
  R->sample( p, bank );
}
