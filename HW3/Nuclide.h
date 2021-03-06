#ifndef _NUCLIDE_HEADER_
#define _NUCLIDE_HEADER_

#include <vector>
#include <string>
#include <memory>

#include "Reaction.h"

class nuclide {
  private:
    std::string nuclide_name;                        // name of nuclide
    double atomicMass;                               // atomic mass fraction A
    std::vector< std::shared_ptr< reaction > > rxn;  // list of reactions
  public:
    nuclide( std::string label ) : nuclide_name(label), atomicMass(0.0) {};   // constructor takes name
    nuclide( std::string label, double A ) : nuclide_name(label), atomicMass(A) {};
    ~nuclide() {};                                   // destructor

    std::string name() { return nuclide_name; }      // return name of nuclide
    std::vector< std::shared_ptr< reaction > > getReactions() {return rxn;} ; // return list of reactions
    void addReaction( std::shared_ptr< reaction > ); // add a reaction to the list of reactions
    double total_xs( particle* p );                  // return the total micro xs
    std::shared_ptr< reaction > sample_reaction( particle* p );   // samples a random reaction based on micro xs
    double getA() { return atomicMass; }             // return atomic mass fraction A
};


#endif
