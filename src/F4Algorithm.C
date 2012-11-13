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
#include "F4Algorithm.H"
#include "Matrix.H"

#include <stdio.h>
#define BREAKPOINT {while(getchar() != '\n');}

using namespace std;

void F4::updatePairs(F4PairSet& pairs, vector<Polynomial>& polys, bool initial) 
{
	double timer = seconds();
	size_t t = groebnerBasis.size();
	for(size_t i = 0; i < polys.size(); i++)
	{
		bool insertIntoG = true;
		Polynomial& h = polys[i];
		// Cancel in P all pairs (i,j) which satisfy T(i,j) = T(i,j,t), T(i,t) != T(i,j) != T(j,t) [ B_t(i,j) ]
		F4PairSet P1(pairs.key_comp());
		for(set<F4Pair>::iterator it = pairs.begin(); it != pairs.end(); it++) 
		{
			if( !it->LCM.isDivisibleBy(h.LT())  || h.lcmLT(groebnerBasis[it->i]) == it->LCM  ||  h.lcmLT(groebnerBasis[it->j]) == it->LCM  ) { 
				P1.insert( *it );
			}   
		}   
		swap(pairs, P1);

		// Let D1 := {(i,t) | 1 <= i < t }.
		for(size_t i = 0; insertIntoG && i < groebnerBasis.size(); i++)
		{
			if(inGroebnerBasis[i] && h.LT().isDivisibleBy(groebnerBasis[i].LT())) 
			{
				insertIntoG = false;
			}
		}

		if(insertIntoG)
		{   
			vector<bool> D1(inGroebnerBasis.begin(), inGroebnerBasis.end());
			// Cancel in D1 each (i,t) for which a (j,t) exists s.t. T(i,t) is a proper multiple of T(j,t) [ M(i,t) ]
			for(size_t i = 0; i < D1.size(); i++) 
			{
				Term a = h.lcmLT(groebnerBasis[i]);
				for(size_t j = 0; D1[i] && j < D1.size(); j++)
				{
					if(i != j && D1[j])				
					{
						Term b = h.lcmLT(groebnerBasis[j]);
						if(a.isDivisibleBy(b) && a != b)
						{
							D1[i] = false;
						}
					}
				}
			}

			// In each nonvoid subset { (j,t) | T(j,t) = tau } ...
			F4PairSet P2(pairs.key_comp());
			for(size_t i = 0; i < D1.size(); i++)
			{
				if(D1[i])
				{
					Term LCM = groebnerBasis[i].lcmLT(h);
					F4Pair newpair( LCM, i, t, LCM == groebnerBasis[i].LT().mul(h.LT()), max(groebnerBasis[i].sugar() - groebnerBasis[i].LT().deg(), h.sugar() - h.LT().deg()) + LCM.deg() );
					pair<set<F4Pair>::iterator,bool> ret;
					// TODO: Efficient ...
					ret = P2.insert( newpair );
					if(newpair.marked && ret.second)
					{
						P2.erase(ret.first);
						P2.insert(newpair);
					}
				}
			}

			// Finally delete all (i,t) with T(i)T(j) = T(i,t).
			for(set<F4Pair>::iterator it = P2.begin(); it != P2.end(); it++)
			{   
				if(!it->marked)
				{ 
					pairs.insert(*it);
				}  
			}

			for(size_t j = 0; j < groebnerBasis.size(); j++)
			{   
				if(inGroebnerBasis[j] && groebnerBasis[j].LT().isDivisibleBy(h.LT()))
				{   
					inGroebnerBasis[j] = false;
				}
			}
		groebnerBasis.push_back( h );
		inGroebnerBasis.push_back( insertIntoG );
		t++;
		}

	}
	//*out << "UPDATE:\t" << seconds() - timer << "\n";
	updateTime += seconds() - timer;
}

struct s
{
	bool operator() (const F4Pair& a, const F4Pair& b)
	{
		return a.sugar < b.sugar;
	}
};


void printPolyMatrix(vector<Polynomial>& v, const TOrdering* O)
{
	Term::comparator t(O, true);
	set<Term, Term::comparator> terms(t);
	for(size_t i = 0; i < v.size(); i++) {
		for(size_t j = 0; j < v[i].size(); j++) {
			terms.insert(v[i][j].second);
		}
	}


	for(size_t i = 0; i < v.size(); i++) {
		set<Term, Term::comparator>::iterator it = terms.begin();
		for(size_t j = 0; j < v[i].size(); j++) {
			for(;*it != v[i][j].second;it++) { std::cout << " 0 "; }
			std::cout << " " << v[i][j].first << " ";
			it++;
		}
		for(;it != terms.end(); it++) { std::cout << " 0 "; }
		std::cout << "\n";
	}

}

PolynomialSequence F4::prepare(F4PairSet& pairs)
{
	//double timer = seconds();
	// SELECTION
	std::cout<<"in F4::prepare, O: "<<O<<std::endl;
	Term::comparator tog(O, true);
	vector<F4Pair> tmp(pairs.begin(), pairs.end());
	sort(tmp.begin(), tmp.end(), s());
	currentDegree = tmp.begin()->sugar;

	PolynomialSequence poly_seq(currentDegree, monoid, tog, field);

	size_t index;
	// Create pivots
	unordered_map<Term, size_t> pivots;

	// selection
	for(index = 0; index < tmp.size() && tmp[index].sugar == currentDegree; index++)
	{
		poly_seq.push_spoly(tmp[index].LCM,
				groebnerBasis[tmp[index].i],
				groebnerBasis[tmp[index].j]);

		//rows.push_back(make_pair(tmp[index].i, tmp[index].LCM));
		//rows.push_back(make_pair(tmp[index].j, tmp[index].LCM));
		//pivots.insert(make_pair(tmp[index].LCM, 2*index));
	}
	pairs.clear();
	pairs.insert(tmp.begin() + index, tmp.end());
	tmp.clear();
	// SELECTION END

	size_t upper = 2*index;

	unordered_set<Term> termsUnordered;

	vector<vector<Monomial> > rightSide;

	//double testtimer = 0;

	// Choose reductor polynomials from those in the current basis
	//for(size_t i = 0; i < rows.size(); i++) 
	for(size_t i = 0; i < poly_seq.size(); i++) 
	{
		//size_t currentRow = rows[i].first;
		//rows[i].second->divToVector(groebnerBasis[currentRow].LT(), ir);
		// For the pivot rows (even rows and lower part) we start at 1 

		// precalculated monomials
		//Term ir = rows[i].second.div(groebnerBasis[currentRow].LT());
		//const Term* ir = rows[i].second->div(groebnerBasis[currentRow].LT());
		//vector<const Term*> pcm = ir->mulAll(groebnerBasis[currentRow], threads, testtimer);
		vector<Term> pcm = poly_seq.get_terms_vector(i);
		vector<coeffType> coeffs = poly_seq.get_coeff_vector(i);

		for(size_t j =  (i > upper || i % 2 == 0 ? 1 : 0);  j < pcm.size() ; j++) 
		{
			//coeffType coeff = coeffs[j];
			const Term t = pcm[j];
			// 50% - inserting in the hash map is expensive
			bool already_considered = termsUnordered.count(t) > 0;	
			
			bool found = false;

			if(!already_considered) {
				// If there is not yet a pivot for t, try to
				// create one
				//found = pivots.count(t) > 0;
				found = poly_seq.check_pivot(t);
				if(!found)
				{
					//for(size_t k = 0; !found && k < groebnerBasis.size(); k++) 
					for(size_t k = 0; k < groebnerBasis.size(); k++) 
					{
						if(inGroebnerBasis[k] && t.isDivisibleBy(groebnerBasis[k].LT())) {
							found = true;
							poly_seq.push_back(t.div(groebnerBasis[k].LT()), groebnerBasis[k], false);
							//rows.push_back(make_pair(k, t));
							//pivots.insert(make_pair(t, rows.size()-1));
						}
					}
				}
			}

			// Eliminate if possible
			if(!found && !already_considered) {
				if(!already_considered) termsUnordered.insert(t);
			}
		}
	}

	//terms.insert(termsUnordered.begin(), termsUnordered.end());
	poly_seq.set_terms<unordered_set<Term> >(termsUnordered.begin(), termsUnordered.end());

	return poly_seq;
}

void F4::reduce(F4PairSet& pairs, vector<Polynomial>& polys)
{
  //double timer = seconds();
  
  
  //vector<vector<size_t> > setOffset; 
  PolynomialSequence poly_seq = prepare(pairs);

  poly_seq.reduce(polys);
   
}

/*
void F4::postReduce(vector<Polynomial>& polys) 
{

}
*/

vector<Polynomial> F4::operator()(vector<Polynomial>& generators, const TOrdering* o, CoeffField* field, int threads, int verbosity, std::ostream& output)
{
  double start = seconds();
  
  this->field = field;
  
	this->threads = threads;
	this->O = o;
	this->verbosity = verbosity;
	this->out = &output;
	F4PairComparator f4pc(o);
	F4PairSet pairs(f4pc);
	updateTime = 0;
	prepareTime = 0;
	reductionTime = 0;

	sort(generators.begin(), generators.end(), PolynomialComparator(o, true));		
	
	//normalize
	for_each(generators.begin(), generators.end(), bind(mem_fn(&Polynomial::normalize), _1, field));

	updatePairs(pairs, generators,true);

	while( !pairs.empty() ) {
		vector<Polynomial> polys;
		reduce(pairs, polys);
		if(!polys.empty()) {
			updatePairs(pairs, polys);
		}
	}
	if(verbosity & 2) {
		*out << "Reduction (s): \t" << reductionTime << "\n";
	}
	if(verbosity & 4) {
		*out << "Prepare (s): \t" << prepareTime << "\n";
	}
	if(verbosity & 8) {
		*out << "Update (s): \t" << updateTime << "\n";
	}
	vector<Polynomial> result;
	for(size_t i = 0; i < groebnerBasis.size(); i++)
	{
		if(inGroebnerBasis[i])
			result.push_back(groebnerBasis[i]);
	}
	if(verbosity & 1) {
		*out << "Runtime (s):\t" << seconds() - start << "\n";
	}
	return result;
}
