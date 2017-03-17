#include "Source.h"

std::stack<particle> source::sample() {
  std::stack<particle> pbank;
  pbank.push( particle( dist_pos->sample(), dist_dir->sample() ) );
  return pbank;
}
