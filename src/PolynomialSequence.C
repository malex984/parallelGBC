#include "../include/PolynomialSequence.H"

void PolynomialSequence::push_back(Polynomial& poly)
{
	seq.push_back(make_pair(monoid.getOne(), poly));
}

void PolynomialSequence::push_back(const Term* t, Polynomial& poly)
{
	seq.push_back(make_pair(t, poly));
}
