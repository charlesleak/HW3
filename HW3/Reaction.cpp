#include "Reaction.h"
#include "Particle.h"
#include "Distribution.h"

void  capture_reaction::sample( particle* p, std::stack< particle >* bank, double A ) {
  // kill the particle and leave the bank unmodified
  p->kill();
}

void  scatter_reaction::sample( particle* p, std::stack< particle >* bank, double A ) {
  // scatter the particle and leave the bank unmodified
  double m = 939.565; // mass of neutron MeV/c2
  double c = 29.98; // speed of light cm/ns
  double muCM = scatter_dist->sample(); // center of mass, needs to be changed to lab frame
  if ( A == 0) { // in this case, A was not set, so don't deal with center of mass stuff
    p->scatter( muCM );
    return;
  }
  double mu0 = ( muCM + A ) / std::sqrt( 1 + A*A + 2 * muCM * A );
  p->scatter( mu0 );
  // v = sqrt(2E/m)
  double vLI = c * std::sqrt( 2.0 * p->getEnergy() / m ); // velocity in the Lab frame Initially
  double vCM = A / ( 1 + A ) * vLI; // velocity in Center of Mass; could be changed for target motion
  double ECM = 0.5 * m * vCM * vCM; // Energy in Center of Mass frame
  double E = ECM + ( p->getEnergy() + 2.0 * muCM * ( A + 1 ) * std::sqrt( p->getEnergy() * ECM ) ) / std::pow( A + 1, 2 ); // final Energy of particle in lab frame
  p->setEnergy( E );
}

void  fission_reaction::sample( particle* p, std::stack< particle >* bank, double A ) {
  // create random number of secondaries from multiplicity distributon and
  // push all but one of them into the bank, and set working particle to the last one
  // if no secondaries, kill the particle

  int n = multiplicity_dist->sample();
  if ( n <= 0 ) {
    p->kill();
  }
  else {
    // bank all but last particle (skips if n = 1)
    for ( int i = 0 ; i < (n - 1) ; i++ ) {
      particle q( p->pos(), isotropic->sample() );
      q.recordCell( p->cellPointer() );
      bank->push( q );
    }
    // set working particle to last one
    particle q( p->pos(), isotropic->sample() );
    q.recordCell( p->cellPointer() );
    *p = q;
  }
}
