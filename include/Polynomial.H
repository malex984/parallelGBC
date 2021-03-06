/*
 *  This file is part of parallelGBC, a parallel groebner basis computation tool.
 *
 *  parallelGBC is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  parallelGBC is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with parallelGBC.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include <iostream>
#include <vector>
#include "../include/Term.H"
#include "../include/CoeffField.H"
#include "../include/TOrdering.H"

using namespace std;

typedef pair<coeffType, const Term*> Monomial;

class Polynomial {
	private:
		std::vector<Monomial> monomials;
		degreeType sugarDegree;

	public:
		// Create empty polynomial (equal to 0)
		Polynomial(degreeType s = 0) : sugarDegree(s) { }

		Polynomial(std::vector<Monomial>& ms) : monomials(ms), sugarDegree(0)
		{ 
			if(ms.size() > 0)
			{ 
				sugarDegree = ms[0].second->deg();
			}
		}

    Polynomial(const Term* t): sugarDegree(t->deg())
    {
			monomials.push_back(make_pair(1, t));
		}

		Polynomial(vector<Monomial>& ms, bool purify) : monomials(ms), sugarDegree(0)
		{ 
			if(ms.size() > 0)
			{ 
				sugarDegree = ms[0].second->deg();
			}
			if(purify)
			{ 
				for(size_t i = 0; i < monomials.size(); i++)
				{ 
					for(size_t j = i+1; j < monomials.size();) {
						if( monomials[i].second == monomials[j].second ) {
							monomials[i].first += monomials[j].first;
							monomials.erase(monomials.begin() + j);
						} else {
							j++; 
						}
					}
				}
			}
		}

		static Polynomial createInstance(const string& s, TMonoid& m, degreeType min = 1);
		static vector<Polynomial> createList(const string& s, TMonoid& m, degreeType min = 1);

		Monomial operator[](size_t i) const {
			return monomials[i];
		}

		size_t size() const {
			return monomials.size();
		}

		size_t capacity() const {
			return monomials.capacity();
		}

		Polynomial mul(const Term* t) const;

		void mulBy(const Term* t);

		bool isZero() const {
			return monomials.size() < 1 || monomials[0].first == 0; // Attention 0*LT+k*T with k != 0 will break this!
		}

		void sub(const Polynomial& other, const TOrdering* O, const CoeffField* f);

		void normalize(const CoeffField* field);

		void mulBy(coeffType l, const CoeffField* f);

		void bringIn(const CoeffField* field, bool normalize = true);

		void order(const TOrdering* O);

		const Term* lcmLT(const Polynomial& other) const {
			return LT()->lcm(other.LT());
		}

		const Term* LT() const {
			return monomials[0].second;
		}

		coeffType LC() const {
			return monomials[0].first;
		}

//		void setSugar(degreeType s) {
//			sugarDegree = s;
//		}

		degreeType sugar() const {
			return sugarDegree;
		}

		void push_back(const Monomial& m) {
			monomials.push_back(m);
		}
};

ostream& operator<< (ostream &out, const Polynomial &poly);
ostream& operator<< (ostream& out, const vector<Polynomial> &polys);

class MonomialComparator {
	public:
		const TOrdering* O;
		greater<coeffType> g;


		MonomialComparator(const TOrdering* O) : O(O) {}

		bool operator() (const Monomial& lhs, const Monomial& rhs) const;
};

class PolynomialComparator {
	public:
		const TOrdering* O;
		bool gt; // if greater

		PolynomialComparator(const TOrdering* O, bool gt = false): O(O), gt(gt) {}

		bool operator() (const Polynomial& lhs, const Polynomial& rhs) const 		{
			return gt ? O->cmp(rhs.LT(), lhs.LT()) < 0 : O->cmp(lhs.LT(), rhs.LT()) < 0;
		}
};



#endif
