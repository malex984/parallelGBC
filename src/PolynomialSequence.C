#include "../include/PolynomialSequence.H"
#include "../include/F4Algorithm.H"


void PolynomialSequence::push_spoly(const Term* lcm, const Polynomial& left,
		const Polynomial& right)
{
	// do not check if the pivot already exists in the list
	// only insert a pivot for one of the parts
	push_back(lcm->div(left.LT()), left, false, false);
	push_back(lcm->div(right.LT()), right, false, true);
	n_spolys++;
}

void PolynomialSequence::push_back(const Polynomial& poly)
{
	push_back(monoid.getOne(), poly, true, true);
}

// utility function to pick the coefficient from a "Monomial"
int firstElement(Monomial &p) {
	return p.first;
}
void PolynomialSequence::push_back(const Term* t, const Polynomial& poly,
		bool check, bool pivot)
{
	// if we check for the existence of a pivot, we are definitely adding
	// a new one
	assert ( !( check && !pivot ) );

	// add polynomial to the sequence
	seq.push_back(make_pair(t, poly));

	double dummy_timer;
	// add the terms of monomial*polynomial
	// we store the terms in a temporary list to make it easier to
	// access the leading term when checking for an already existing pivot
	const terms_type cterms = t->mulAll(poly, F4_THREADS, dummy_timer );
	seq_terms.push_back(cterms);

	coeffs_type coeffs(poly.size());
	std::transform(poly.begin(), poly.end(), std::back_inserter(coeffs),
			firstElement);
	seq_coeffs.push_back(coeffs);
	if (pivot) {
		size_t pivot_ind = seq.size() - 1;
		pivots.insert(make_pair(t->mul(poly.LT()), pivot_ind));
		if (check) {
			pivots_type::iterator prev_pivot_itr =
				pivots.find(cterms[0]);
			if (prev_pivot_itr != pivots.end())
			{
				// TODO: select the short one as the pivot
				pivot_ops[t].push_back(
						make_pair(pivot_ind, coeffs[0]) );
			}
		}
	}
}

const std::vector<coeffType>& PolynomialSequence::get_coeff_vector(size_t ind) const 
{
	return seq_coeffs[ind];
}

const std::vector<const Term*>& PolynomialSequence::get_terms_vector(size_t ind) const
{
	return seq_terms[ind];
}

void PolynomialSequence::reduce(vector<Polynomial>& polys) const
{
	// terms must be set before calling reduce
	assert(terms.size() > 0);
	size_t nrows = seq.size();
	Matrix pmatrix(nrows, terms.size());
	initialize_matrix(pmatrix);
	// arrange eliminations in pivot columns to be parallel friendly
	vector<vector<F4Operation> > ops;
	order_pivot_ops(ops);


	// apply pivot_ops to the matrix
	pReduce(pmatrix, ops);

	//cerr << "pReduced matrix: " << *rs << endl;

	// this will be used to mark zero rows in the row reduced form
	vector<bool> empty_rows(2*n_spolys, false); // too large, FIX?

	gauss(pmatrix, empty_rows);

	//cerr << "Gaussed matrix: " << *rs << endl;  

	get_new_polys(polys, pmatrix, empty_rows);
}

void PolynomialSequence::get_new_polys(vector<Polynomial> &polys, Matrix& pmatrix, vector<bool> empty_rows) const
{
	// 	for(size_t i = 1; i < upper; i+=2)
	for (size_t i = 1; i < 2*n_spolys; i += 2)
	{
		if(!empty_rows[i])
		{
			Polynomial p(sugar_degree);
			size_t j = 0;
			for(terms_type::const_iterator it = terms.begin();
					it != terms.end(); it++, j++) 
			{
				const coeffType v = pmatrix.getEntry(i, j);

				if(v != 0)
					p.push_back(make_pair(v, *it));

			}
			polys.push_back( p );
		}
	}

	//cerr << "Resulting polynomials: " << polys.size() << " new elements: \n" << polys << endl;
}

void PolynomialSequence::initialize_matrix (Matrix& pmatrix) const
{
	//   size_t pad = __COEFF_FIELD_INTVECSIZE;


	//Matrix prs = Matrix::allocate(nrows, (( terms.size()+pad-1 )/ pad ) * pad, 0); 
	seq_coeffs_type::const_iterator cur_coeffs = seq_coeffs.begin();
	size_t i = 0;
	for(seq_terms_type::const_iterator cur_terms = seq_terms.begin();
			cur_coeffs != seq_coeffs.end();
			cur_coeffs++, cur_terms++, i++)
	{
		terms_type::const_iterator cur_term = (*cur_terms).begin();
		coeffs_type::const_iterator cur_coeff = (*cur_coeffs).begin();
		size_t j = 0;
		for(terms_type::const_iterator it = terms.begin();
				cur_coeff != (*cur_coeffs).end() ; it++, j++)
		{
			if(*cur_term == *it)
			{
				pmatrix.setEntry(i, j, *cur_coeff);
			}
		}
	}
}

void PolynomialSequence::order_pivot_ops(vector<vector<F4Operation> > &ops) const
{
	map<const Term*, vector<pair<size_t, coeffType> >, TermComparator>
		pivotOpsOrdered(pivot_ops.begin(), pivot_ops.end(),
				term_comparator);
	ops.push_back( vector<F4Operation >() );

	#if 1 
	vector<size_t> l(seq.size(), 0);
	for(map<const Term*, vector<pair<size_t, coeffType> >, TermComparator>::reverse_iterator it = pivotOpsOrdered.rbegin(); it != pivotOpsOrdered.rend(); it++)
	{
		size_t o = pivots.find(it->first)->second;
		for(size_t i = 0; i < it->second.size(); i++)
		{
			size_t t = it->second[i].first;
			if(l[ o ] > l[t]) {
				l[t] = l[o];
			}
			ops[ l[t] ].push_back( F4Operation(t,o,it->second[i].second) );
			l[t]++; // one operation per level, attention this also affects the following if-statements

			if(l[t] >= ops.size()) {
				ops.push_back( vector<F4Operation>() );
			}
		}
	}
	#else
	size_t l = 0;
	for(map<const Term*, vector<pair<size_t, coeffType> >, TermComparator>::reverse_iterator it = pivotOpsOrdered.rbegin(); it != pivotOpsOrdered.rend(); it++)
	{
		size_t o = pivots[it->first];
		for(size_t i = 0; i < it->second.size(); i++)
		{
			size_t t = it->second[i].first;
			ops[ l ].push_back( F4Operation(t, o, it->second[i].second) );
		}
		l++;
		ops.push_back( vector<F4Operation>() );
	}
	#endif
	pivotOpsOrdered.clear();
	ops.pop_back();
	//prepareTime += seconds() - timer;
	//cout << "Operations:\t" << ops.size() << "\n";
	//cout << "Preparation:\t" << seconds() - timer << "\n";
}

void PolynomialSequence::gauss(Matrix& pmatrix, vector<bool>& empty) const
{
	pmatrix.gauss(2*n_spolys, empty, field);
}

void PolynomialSequence::pReduce(Matrix& pmatrix,
		vector<vector<F4Operation> >& ops) const
{
	for(size_t i = 0; i < ops.size(); i++)
	{
		//		#pragma omp parallel for num_threads( threads ) 
		for(size_t j = 0; j < ops[i].size(); j++)
			pmatrix.pReduce(ops[i][j].target, ops[i][j].oper,
					ops[i][j].factor, field);
	}
}

