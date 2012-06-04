#include "../include/Term.H"
#include "../include/Polynomial.H"

#include <vector>
#include <utility>

class PolynomialSequence {
public:
	PolynomialSequence(const TMonoid& monoid) : monoid(monoid),
		seq(), terms() {};
	~PolynomialSequence();

	void push_back(Polynomial& poly);
	void push_back(const Term* t, Polynomial& poly);
	template<class container>
	void set_terms(typename container::const_iterator beg,
			typename container::const_iterator end) {
		terms.insert(terms.begin(), beg, end);
	};

private:
	const TMonoid& monoid;
	std::vector<std::pair<const Term*, Polynomial> > seq;
	std::vector<const Term*> terms;
};

