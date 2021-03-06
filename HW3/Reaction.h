#ifndef _REACTION_HEADER_
#define _REACTION_HEADER_

#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include <stack>
#include <utility>
#include <cmath>

#include "Particle.h"
#include "Distribution.h"

class reaction {
  private:
    double rxn_xs;
    double oneOverVFactor; // rxn_xs + k * 1/sqrt(E)
  protected:
    std::string rxn_name;
  public:
    reaction( double x, double k ) : rxn_xs(x), oneOverVFactor(k) {};
    ~reaction() {};

    virtual std::string name() final { return rxn_name; };
    virtual double xs( particle* p ) final { 
      return rxn_xs + oneOverVFactor / std::sqrt( p->getEnergy() );
    };
    virtual void sample( particle* p, std::stack<particle>* bank, double A ) = 0;   // pure virtual
};

class capture_reaction : public reaction {
  private:
 
  public:
    capture_reaction( double x, double k ) : reaction(x,k) { rxn_name = "capture"; }; //construct with xs
    ~capture_reaction() {};

    void sample( particle* p, std::stack<particle>* bank, double A );     // sample capture
};

class scatter_reaction : public reaction {
  private:
    std::shared_ptr< distribution<double> > scatter_dist; 
  public:
    scatter_reaction( double x, double k, std::shared_ptr< distribution<double> > D ) : // construct with xs and angular distribution
      reaction(x,k), scatter_dist(D) { rxn_name = "scatter"; };
    ~scatter_reaction() {};

    void sample( particle* p, std::stack<particle>* bank, double A );         // sample scatter, including Energy change
};

class fission_reaction : public reaction {
  private:
    std::shared_ptr< distribution<int> >   multiplicity_dist; 
    std::shared_ptr< distribution<point> > isotropic;
  public:
    fission_reaction( double x, double k, std::shared_ptr< distribution<int> > D ) : // construct with xs and multiplicity distribution
    reaction(x,k), multiplicity_dist(D) { 
         rxn_name = "fission";
         isotropic = std::make_shared< isotropicDirection_distribution > ( "isotropic" ); 
       };
    ~fission_reaction() {};

    void sample( particle* p, std::stack<particle>* bank, double A );      // sample fission
};

#endif
