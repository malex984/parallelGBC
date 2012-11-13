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
#ifndef F4_ALGORITHM_H
#define F4_ALGORITHM_H

#include "Term.H"
#include "Polynomial.H"
#include "PolynomialSequence.H"
#include "F4Utils.H"

#include <set>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <boost/bind.hpp>
#include <iostream>

#ifdef F4_THREADS
	#define __F4_ALGORITHM_THREADS F4_THREADS
#else
	#define __F4_ALGORITHM_THREADS 1
#endif

class F4Pair {
	public:
		Term LCM;
		size_t i;
		size_t j;
		bool marked;
		degreeType sugar;

		F4Pair(const Term LCM, size_t i, size_t j, bool marked, degreeType sugar) : LCM(LCM), i(i), j(j), marked(marked), sugar(sugar)	 {};
};

class F4PairComparator {
	public:
		const TOrdering* O;

		F4PairComparator(const TOrdering* O) : O(O) {}

		bool operator() (const F4Pair& lhs, const F4Pair& rhs) const 
		{
			return O->cmp(lhs.LCM, rhs.LCM) < 0;
		}
};

typedef std::set<F4Pair,F4PairComparator> F4PairSet;

class Matrix;

class F4 {
	public:
    CoeffField* field;
    TMonoid& monoid;

		std::vector<Polynomial> groebnerBasis;
		std::vector<bool> inGroebnerBasis;
		const TOrdering* O;
		degreeType currentDegree;
		int threads;
		/* Verbosity, which can be changed during runtime, nothing which should
		 * influence performance.
		 *
		 * 1 - Runtime
		 * 2 - Reduction time
		 * 4 - Prepare time
		 * 8 - Update time
		 * 16 - Print sugar degree during reduction step
		 * 32 - Print time of reduction step
		 * 64 - Print matrix size during reduction step 
		 */
		int verbosity;
		// There we write verbose messages
		std::ostream* out;

		F4(TMonoid& monoid) : monoid(monoid) {}
		

		// statistics
		double reductionTime;
		double prepareTime;
		double updateTime;

    // select() is part of reduce
		
		//void gauss(Matrix* pm, size_t upper, vector<bool>& empty);
		//void pReduce(vector<vector<F4Operation> >& ops, Matrix* prs);
		PolynomialSequence prepare(F4PairSet& pairs);
		void reduce(F4PairSet& pairs, std::vector<Polynomial>& polys);
		//size_t prepare(F4PairSet& pairs, vector<Polynomial>& polys, vector<vector<F4Operation> >& ops, set<const Term*, TermComparator>& terms, Matrix **prs);
		void postReduce(std::vector<Polynomial>& polys);
		void updatePairs(F4PairSet& pairs, std::vector<Polynomial>& polys, bool initial = false);
		std::vector<Polynomial> operator()(std::vector<Polynomial>& generators, const TOrdering* O, CoeffField* field, int threads = __F4_ALGORITHM_THREADS, int verbosity = 0, std::ostream& output = std::cout);
};

#endif
