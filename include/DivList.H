#include "mathic/mathic.h"
#include "../include/Term.H"

#include <iostream>

class BasicConfiguration {
public: // the methods that Mathic calls must be public

  // The configuration needs to know how many variables there are
  // so that it can tell Mathic. Variables are indexed from 0 to
  // numberOfVariables - 1.
  BasicConfiguration(size_t numberOfVariables):
    varCount(numberOfVariables) {}

  // The type of an exponent in a monomial.
  typedef degreeType Exponent;

  // The type of monomials.
  typedef const Term* Monomial;

  // The type of monomials that are entries in the data structure.
  // For some applications we might want this to be a completely
  // different type from Monomial, but for this simple example we will
  // just let this be the same type as Monomial.
  typedef std::pair<Monomial, size_t> Entry;

  // This method is how Mathic determines how many variables there are.
  // Mathic does not store this number and instead calls this method every time
  // it needs to know this information, so this method should run quickly.
  size_t getVarCount() const {
    return varCount;
  }

  // Mathic uses this function to extract an Exponent from a Monomial.
  Exponent getExponent(Monomial m, size_t var) const {
    MATHIC_ASSERT(var < m.size());
    return (*m)[var];
  }
  Exponent getExponent(Entry m, size_t var) const {
    return (*m.first)[var];
  }


  // Mathic uses this function to determine if a divides b.
  // This method can be implemented using the other methods on a
  // configuration as done here. The point of this method is that
  // some representations of monomials allow faster divisiblity
  // checks than is achieved by extracting every exponent.
  //
  // We only need one method called divides since Monomial and Entry
  // are the same type. Otherwise we would need 3 versions
  // taking parameters (Monomial, Entry), (Entry, Monomial)
  // and (Monomial,Monomial). Thankfully here we only need
  // (Monomial,Monomial).
  bool divides(const Entry& a, const Monomial& b) const {
	  return b->isDivisibleBy(a.first);
  }
  bool divides(const Monomial& a, const Entry& b) const {
	  return b.first->isDivisibleBy(a);
  }
  bool divides(const Monomial& a, const Monomial& b) const {
	  return b->isDivisibleBy(a);
  }

private:
  size_t varCount; // the number of variables
};

// We can now combine the BasicConfiguration with options for a KDTree
// to get a configuration type for a KDTree.
class KDTreeConfiguration :
  public BasicConfiguration, // sets std::vector as the monomial representation
  public mathic::KDTreeSuggestedOptions // sets default options for the KDTree
{
public:
  KDTreeConfiguration(size_t varCount): BasicConfiguration(varCount) {}
};

// With a configuration type, we can now create a KDTree type that uses
// the std::vector representation and using default options.
class MyKDTree : public mathic::KDTree<KDTreeConfiguration> {
public:
  MyKDTree(size_t varCount): 
    KDTree<KDTreeConfiguration>(KDTreeConfiguration(varCount)) {}
};
