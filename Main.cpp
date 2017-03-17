#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <cmath>
#include <vector>
#include <iomanip>
#include <cfloat>
#include <limits>
#include <memory>
#include <cassert>
#include <time.h>
#include <stack>

#include "pugixml.hpp"
#include "Distribution.h"
#include "Particle.h"
#include "Surface.h"
#include "Cell.h"
#include "Source.h"
#include "Reaction.h"
#include "Material.h"
#include "Nuclide.h"
#include "Point.h"

// NOTE: Transport starts at line 425

template< typename T >
std::shared_ptr< T > findByName( std::vector< std::shared_ptr< T > > vec, std::string name ) {
  for ( auto v : vec ) {
    if ( v->name() == name ) { return v; }
  }
  return nullptr;
}

// Find what cell a particle is in
void FindCurrentCell ( particle* p , std::vector < std::shared_ptr < cell > >* vec ) {
    std::shared_ptr< cell > cel = nullptr;
    for ( int v = 0; v < vec->size(); v++ ) { // std::shared_ptr<cell> v : vec
        if ( vec->at(v)->testPoint( p->pos() ) ) { // v->testPoint(p->pos())
            cel = vec->at(v);
            break;
        } 
    }
    assert( cel ); // terminates program if cell is null
    p->recordCell( cel );
} // end find cell function


int main() 
{
    //*
    int long long NSamples = 1e7;

    double err = 1000.0 * std::numeric_limits < double > ::epsilon();

//*
    // user enters the XML file name and pugixml will attempt to load
  std::string input_file_name;
  std::cout << " Enter XML input file name: " << std::endl;
  std::cin  >> input_file_name;

  pugi::xml_document input_file;
  pugi::xml_parse_result load_result = input_file.load_file( input_file_name.c_str() );

  // check to see if result failed and throw an exception if it did
  if ( ! load_result ) {
    std::cout << load_result.description() << std::endl;
    throw;
  }

  // distributuions
  std::vector< std::shared_ptr< distribution<double> > >  double_distributions;
  std::vector< std::shared_ptr< distribution<int>    > >  int_distributions;
  std::vector< std::shared_ptr< distribution<point>  > >  point_distributions;
  pugi::xml_node input_distributions = input_file.child("distributions");

  // find total number of distributions
  int num_distributions = 0;
  for ( auto d : input_distributions ) { num_distributions++; }

  // since distributions may depend on other distributions, need to iterate
  int set_distributions = 0;
  while ( set_distributions < num_distributions ) {
    int previous_set_distributions = set_distributions;

    for ( auto d : input_distributions ) {
      std::string type = d.name();
      std::string name = d.attribute("name").value();
      std::string data = d.attribute("datatype").value();

      if ( data == "double" ) {
        // skip rest of loop if distribution already done
        if ( findByName( double_distributions, name ) ) { continue; }

        std::shared_ptr< distribution<double> > Dist;
        if ( type == "delta" ) {
          double a = d.attribute("a").as_double();
          Dist = std::make_shared< arbitraryDelta_distribution< double > > ( name, a );
        }
        else if ( type == "uniform" ) {
          double a = d.attribute("a").as_double();
          double b = d.attribute("b").as_double();
          Dist = std::make_shared< uniform_distribution > ( name, a, b );
        }
        else if ( type == "linear" ) {
          double a  = d.attribute("a").as_double();
          double b  = d.attribute("b").as_double();
          double fa = d.attribute("fa").as_double();
          double fb = d.attribute("fb").as_double();
          Dist = std::make_shared< linear_distribution > ( name, a, b, fa, fb );
        }
        else if ( type == "henyeyGreenstein" ) {
          double a = d.attribute("a").as_double();
          Dist = std::make_shared< HenyeyGreenstein_distribution > ( name, a );
        }
        else {
          std::cout << "unsupported distribution with data type " << data << std::endl;
          throw;
        }
        double_distributions.push_back( Dist );
      }
      // integer-valued distributions
      else if ( data == "int" ) {
        // skip rest of loop if distribution already done
        if ( findByName( int_distributions, name ) ) { continue; }

        std::shared_ptr< distribution<int> > Dist;
        if ( type == "delta" ) {
          double a = d.attribute("a").as_int();
          Dist = std::make_shared< arbitraryDelta_distribution< int > > ( name, a );
        }
        else if ( type == "meanMultiplicity" ) {
          double nubar = d.attribute("nubar").as_double();
          Dist = std::make_shared< meanMultiplicity_distribution > ( name, nubar );
        }
        else if ( type == "terrellFission" ) {
          double nubar = d.attribute("nubar").as_double();
          double sigma = d.attribute("sigma").as_double();
          double b     = d.attribute("b").as_double();
          Dist = std::make_shared< TerrellFission_distribution > ( name, nubar, sigma, b );
        }
        else {
          std::cout << "unsupported distribution with data type " << data << std::endl;
          throw;
        }
        int_distributions.push_back( Dist );
      }
      else if ( data == "point" ) {
        // skip rest of loop if distribution already done
        if ( findByName( point_distributions, name ) ) { continue; }

        std::shared_ptr< distribution< point > > Dist;
        if ( type == "delta" ) {
          double x = d.attribute("x").as_double(); 
          double y = d.attribute("y").as_double(); 
          double z = d.attribute("z").as_double();         
          Dist = std::make_shared< arbitraryDelta_distribution< point > > ( name, point( x, y, z ) );
        }
        else if ( type == "isotropic" ) {
          Dist = std::make_shared< isotropicDirection_distribution > ( name );
        }
        else if ( type == "anisotropic" ) {
          double u = d.attribute("u").as_double(); 
          double v = d.attribute("v").as_double(); 
          double w = d.attribute("w").as_double();         
          std::shared_ptr< distribution<double> > angDist = 
            findByName( double_distributions, d.attribute("distribution").value() );
      
          // in the angular distribution does not yet, skip to the end of the loop
          if ( ! angDist ) { continue; }

          Dist = std::make_shared< anisotropicDirection_distribution > ( name, point( u, v, w ), angDist );
        }
        else if ( type == "independentXYZ" ) {
          std::shared_ptr< distribution<double> > distX = findByName( double_distributions, d.attribute("x").value() ); 
          std::shared_ptr< distribution<double> > distY = findByName( double_distributions, d.attribute("y").value() ); 
          std::shared_ptr< distribution<double> > distZ = findByName( double_distributions, d.attribute("z").value() ); 

          // if any of these distributions have not yet been resolved, skip to the end of the loop
          if ( !distX || !distY || !distZ ) { continue; }

          Dist = std::make_shared< independentXYZ_distribution > ( name, distX, distY, distZ );
        }
        else {
          std::cout << "unsupported " << data << " distribution of type " << type << std::endl;
          throw;
        }
        point_distributions.push_back( Dist );
      }
      else {
        std::cout << "unsupported distribution with data type " << data << std::endl;
        throw;
      }
      // if we reach here, assume distribution has been set
      set_distributions++;
    }
    // check to see if number of distributions has increased, if not, caught in an infinite loop
    if ( previous_set_distributions == set_distributions ) { 
      std::cout << "distributions could not be resolved. " << std::endl;
      throw;
    }
  }

  // iterate over nuclides
  std::vector< std::shared_ptr<nuclide> > nuclides;
  pugi::xml_node input_nuclides = input_file.child("nuclides");
  for ( auto n : input_nuclides ) {
    std::string name = n.attribute("name").value();

    std::shared_ptr< nuclide > Nuc = std::make_shared< nuclide > ( n.attribute("name").value() );
    nuclides.push_back( Nuc );

    // iterate over its reactions
    for ( auto r : n.children() ) {
      std::shared_ptr< reaction > Rxn;
      std::string rxn_type = r.name();

      double xs = r.attribute("xs").as_double();
      if ( rxn_type == "capture" ) {
        Nuc->addReaction( std::make_shared< capture_reaction > ( xs ) );
      }
      else if ( rxn_type == "scatter" ) {
        std::string dist_name = r.attribute("distribution").value();
        std::shared_ptr< distribution<double> > scatterDist = findByName( double_distributions, dist_name );
        if ( scatterDist ) {
          Nuc->addReaction( std::make_shared< scatter_reaction > ( xs, scatterDist ) );
        }
        else {
          std::cout << " unknown scattering distribution " << dist_name << " in nuclide " << name << std::endl;
          throw;
        }
      }
      else if ( rxn_type == "fission" ) {
        std::string mult_dist_name = r.attribute("multiplicity").value();
        std::shared_ptr< distribution<int> > multDist = findByName( int_distributions, mult_dist_name );
        if ( multDist ) {
          Nuc->addReaction( std::make_shared< fission_reaction > ( xs, multDist ) );
        }
        else {
          std::cout << " unknown multiplicity distribution " << mult_dist_name << " in nuclide " << name << std::endl;
          throw;
        }
      }
      else {
        std::cout << "unknown reaction type " << rxn_type << std::endl;
        throw;
      }
    }
  } 

  // iterate over materials
  std::vector< std::shared_ptr<material> > materials;
  pugi::xml_node input_materials = input_file.child("materials");
  for ( auto m : input_materials ) {
    std::string name = m.attribute("name").value();
    double      aden = m.attribute("density").as_double();
    
    std::shared_ptr< material > Mat = std::make_shared< material > ( name, aden );    
    materials.push_back( Mat );

    // iterate over nuclides
    for ( auto n : m.children() ) {
      if ( (std::string) n.name() == "nuclide" ) {
        std::string nuclide_name = n.attribute("name").value();
        double      frac         = n.attribute("frac").as_double();
        
        Mat->addNuclide( findByName( nuclides, nuclide_name ), frac );
      }
    }
  }

  // iterate over surfaces
  std::vector< std::shared_ptr< surface > > surfaces;
  pugi::xml_node input_surfaces = input_file.child("surfaces");
  for ( auto s : input_surfaces ) {
    std::string type = s.name();

    std::shared_ptr< surface > S;
    if ( type == "plane" ) {
      std::string name = s.attribute("name").value();
      double      a    = s.attribute("a").as_double();
      double      b    = s.attribute("b").as_double();
      double      c    = s.attribute("c").as_double();
      double      d    = s.attribute("d").as_double();
      S = std::make_shared< plane > ( name, a, b, c, d );
    }
    else {
      std::cout << " unkown surface type " << type << std::endl;
      throw;
    }

    if ( (std::string) s.attribute("bc").value() == "reflect" ) {
      S->makeReflecting();
    }
    surfaces.push_back( S );
  }

  // iterate over cells
  std::vector< std::shared_ptr< cell > > cells;
  pugi::xml_node input_cells = input_file.child("cells");
  for ( auto c : input_cells ) {
    std::string name = c.attribute("name").value();

    std::shared_ptr< cell > Cel = std::make_shared< cell > ( name );
    cells.push_back( Cel );

    // cell material
    if ( c.attribute("material") ) {
      std::shared_ptr< material > matPtr = findByName( materials, c.attribute("material").value() );
      if ( matPtr ) {
        Cel->setMaterial( matPtr );
      }
      else {
        std::cout << " unknown material " << c.attribute("material").value() << " in cell " << name << std::endl;
        throw;
      } 
   }

    // cell importance
    if ( c.attribute("importance") ) {
      Cel->setImportance( c.attribute("importance").as_double() );
    }
   
    // iterate over surfaces
    for ( auto s : c.children() ) {
      if ( (std::string) s.name() == "surface" ) {
        std::string name  = s.attribute("name").value();
        int         sense = s.attribute("sense").as_int();

        std::shared_ptr< surface > SurfPtr = findByName( surfaces, name );
        if ( SurfPtr ) {
          Cel->addSurface( findByName( surfaces, name ), sense );
        }
        else {
          std::cout << " unknown surface with name " << name << std::endl;
          throw;
        }
      }
      else {
        std::cout << " unknown data type " << s.name() << " in cell " << name << std::endl;
        throw;
      }
    } 
  }

  // iterate over estimatators
  std::vector< std::shared_ptr< estimator > > estimators;
  pugi::xml_node input_estimators = input_file.child("estimators");
  for ( auto e : input_estimators ) {
    std::string type = e.name();
    std::string name = e.attribute("name").value();
    
    std::shared_ptr< estimator > Est;
    if ( type == "current" ) {
      Est = std::make_shared< surface_current_estimator > ( name );

      // get the surfaces
      for ( auto s : e.children() ) {
        if ( (std::string) s.name() == "surface" ) {
          std::string name = s.attribute("name").value();
          std::shared_ptr< surface > SurfPtr = findByName( surfaces, name );
          if ( SurfPtr ) {
            SurfPtr->attachEstimator( Est );
          }
          else {
            std::cout << " unknown surface label " << name << " in estimator " << e.attribute("name").value() << std::endl;
          }
        }
      } 
    }
    else if ( type == "countingSurface" ) {
      Est = std::make_shared< counting_estimator > ( name );

      // get the surfaces
      for ( auto s : e.children() ) {
        if ( (std::string) s.name() == "surface" ) {
          std::string name = s.attribute("name").value();
          std::shared_ptr< surface > SurfPtr = findByName( surfaces, name );
          if ( SurfPtr ) {
            SurfPtr->attachEstimator( Est );
          }
          else {
            std::cout << " unknown surface label " << name << " in estimator " << e.attribute("name").value() << std::endl;
          }
        }
      } 
    }
    else {
      std::cout << "unknown estimator type " << name << std::endl;
      throw;
    }
    
    estimators.push_back( Est );
  }

  // create source
  pugi::xml_node input_source = input_file.child("source");
  pugi::xml_node input_source_position  = input_source.child("position");
  pugi::xml_node input_source_direction = input_source.child("direction");

  std::string pos_dist_name = input_source_position.attribute("distribution").value();
  std::string dir_dist_name = input_source_direction.attribute("distribution").value();

  std::shared_ptr< distribution< point > > posDist = findByName( point_distributions, pos_dist_name );
  std::shared_ptr< distribution< point > > dirDist = findByName( point_distributions, dir_dist_name );

  std::shared_ptr< source > src;  
  if ( posDist && dirDist ) {
    src = std::make_shared< source > ( posDist, dirDist );  
  }
  else {
    if ( ! posDist ) { std::cout << " unknown position distribution "  << pos_dist_name << " in source " << std::endl; }
    if ( ! dirDist ) { std::cout << " unknown direction distribution " << dir_dist_name << " in source " << std::endl; }
    throw;
  }

  //*/






    clock_t init, final;
    init=clock();
    int count = 0;

    // Begin loop over histories
    for ( int n = 1; n <= NSamples; n++ ) {

        std::stack < particle > bank = src->sample();
        // std::stack < particle > bank = initSource->sample();


        while ( ! bank.empty() ) {
            particle p = bank.top() ; bank.pop();
            // find cell particle is in
            FindCurrentCell ( &p , &cells );
            std::shared_ptr < cell > currentCell = p.cellPointer();
            // std::cout << currentCell->name() << std::endl;


            while ( p.alive() ) {

                // find distance to nearest boundary
                std::pair < std::shared_ptr < surface > , double > rayIntersect = 
                    currentCell->surfaceIntersect( p.getRay() );
                double distToBound  = rayIntersect.second;
                // std::cout << distToBound << std::endl;

                // find distance to collision
                double distToCollision = -std::log( Urand() )/(currentCell->macro_xs());
                // std::cout << distToCollision << std::endl;

                double transDist = std::fmin( distToCollision , distToBound );
                // move particle to new location
                currentCell->moveParticle( &p , transDist );

                // determine if boundary interaction
                if ( transDist == distToBound ) { // boundary interaction
                    // advance to surface and apply crossSurface method from Surface class
                    rayIntersect.first->crossSurface( &p );
                    
                    FindCurrentCell ( &p , &cells );
                      currentCell = p.cellPointer();
                    // kill particle if region has importace of zero
                    if ( currentCell->getImportance() <= 0.0 ) {
                        p.kill();
                    }


                }
                else { // collision
                    // sample collision
                    currentCell->sampleCollision( &p , &bank );

                    // rayIntersect.first->scoreEstimators( &p );

                }

            } // end while particle alive loop

            
            // current->endHistory();

        } // end while bank !empty loop
        for (auto e: estimators) { e->endHistory();count++;}




    } // end for loop over number of histories

    for (auto e: estimators) { e->report();}
    // current->report();

    final=clock()-init;
    std::cout << "Run time: " << (double)final / ((double)CLOCKS_PER_SEC) << " seconds" << std::endl;
    std::cout << count << std::endl;

    return 0;
}


