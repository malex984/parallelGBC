#include "PolynomialSequence.H"
#include "F4Algorithm.H"

#include <functional>

using std::bind;
using std::string;
using std::function;

void PolynomialSequence::push_spoly(const Term& lcm, const Polynomial& left,
		const Polynomial& right)
{
	// do not check if the pivot already exists in the list
	// only insert a pivot for one of the parts
	push_back(lcm.div(left.LT()), left, false, true);
	push_back(lcm.div(right.LT()), right, true, true);
	n_spolys++;
}

void PolynomialSequence::push_back(const Polynomial& poly)
{
	std::cout<<"Pushing poly "<<poly<<std::endl;
	seq.push_back(std::make_pair(monoid.one_term(), poly));

	push_terms(poly);
	push_coeffs(poly);

	add_pivot(poly.LT(), poly.LC(), true);
}

void PolynomialSequence::push_back(const Term& t, const Polynomial& poly,
		bool check, bool pivot)
{
	std::cout<<"Pushing poly "<<t<<" * "<<poly<<std::endl;
	// if we check for the existence of a pivot, we are definitely adding
	// a new one
	assert ( !( check && !pivot ) );

	// add polynomial to the sequence
	seq.push_back(std::make_pair(t, poly));

	push_terms(t, poly);
	push_coeffs(poly);

	if (pivot) {
		add_pivot(t.mul(poly.LT()), poly.LC(), check);
	}
}
void PolynomialSequence::push_coeffs(const Polynomial& poly) {
	coeffs_type coeffs;
	coeffs.insert(coeffs.begin(), poly.begin_coeffs(), poly.end_coeffs());
	seq_coeffs.push_back(coeffs);
}

// utility function to multiply two terms
// this is used in the lambda expression in push_terms(t, poly)
Term  mul_terms(const Term& t, const Term& u) {
	return t.mul(u);
}

void PolynomialSequence::push_terms(const Term& t, const Polynomial& poly) {
	// add the terms of monomial*polynomial
	// we store the terms in a temporary list to make it easier to
	// access the leading term when checking for an already existing pivot
	////terms_type cterms(poly.size());
	
	// giving a size argument to the constructor requires a default
	// constructor for Term
	terms_type cterms;
	
	// mul_by_t(x) = mul_terms(x, t);
	function<Term(Term)> mul_by_t = bind(mul_terms, _1, t);

	std::transform(poly.begin_terms(), poly.end_terms(),
			std::back_inserter(cterms),
			mul_by_t);
	std::cout<<"Pushing terms with size: "<<cterms.size()<<std::endl;
	seq_terms.push_back(cterms);

}

void PolynomialSequence::push_terms(const Polynomial& poly) {
	// giving a size argument to the constructor requires a default
	// constructor for Term
	terms_type cterms; 
	cterms.insert(cterms.begin(), poly.begin_terms(), poly.end_terms());
	std::cout<<"Pushing terms with size: "<<cterms.size()<<std::endl;
	seq_terms.push_back(cterms);
}

void PolynomialSequence::add_pivot(const Term& t, coeffType coeff, bool check) {
	std::cout<<"Adding pivot "<<t<<std::endl;
	size_t pivot_ind = seq.size() - 1;
	// FIXME: what is this supposed to do?
	// shouldn't it remove previous pivots?
	// what happens to the operations recorded for those pivots?
	if (check) {
		pivots_type::iterator prev_pivot_itr =
			pivots.find(t);
		if (prev_pivot_itr != pivots.end())
		{
			// TODO: select the short one as the pivot
			std::cout<<"pivot op: ind: "<<pivot_ind<<" coeff: "<<coeff<<" prev_ind: "<<prev_pivot_itr->second<<std::endl;
			pivot_ops[t].push_back(
					std::make_pair(pivot_ind, coeff) );
		}
	}
	pivots.insert(std::make_pair(t, pivot_ind));
}

const std::vector<coeffType>& PolynomialSequence::get_coeff_vector(size_t ind) const 
{
	return seq_coeffs[ind];
}

const std::vector<Term>& PolynomialSequence::get_terms_vector(size_t ind) const
{
	return seq_terms[ind];
}

void PolynomialSequence::reduce(std::vector<Polynomial>& polys) const
{
	// terms must be set before calling reduce
	assert(terms.size() > 0);
	size_t nrows = seq.size();
	Matrix pmatrix(nrows, terms.size());
   
	initialize_matrix(pmatrix);
	std::cerr << "matrix: " << pmatrix << std::endl;
	// arrange eliminations in pivot columns to be parallel friendly
	std::vector<F4Operations> ops;
	std::vector<size_t> deps;
	order_pivot_ops(ops, deps);


	// apply pivot_ops to the matrix
	pReduce(pmatrix, ops);

	std::cerr << "pReduced matrix: " << pmatrix << std::endl;

	// this will be used to mark zero rows in the row reduced form
	std::vector<bool> empty_rows(2*n_spolys, false); // too large, FIX?

	gauss(pmatrix, empty_rows);

	std::cerr << "Gaussed matrix: " << pmatrix << std::endl;

	get_new_polys(polys, pmatrix, empty_rows);
}

void PolynomialSequence::get_new_polys(std::vector<Polynomial> &polys,
		Matrix& pmatrix, std::vector<bool> empty_rows) const
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
					p.push_back(std::make_pair(v, *it));

			}
			polys.push_back( p );
		}
	}

	//std::cerr << "Resulting polynomials: " << polys.size() << " new elements: \n" << polys << std::endl;
}

void PolynomialSequence::initialize_matrix (Matrix& pmatrix) const
{
	//   size_t pad = __COEFF_FIELD_INTVECSIZE;

//         std::cerr << "terms: " << terms << std::endl;


	assert(seq_terms.size() == seq_coeffs.size());
	//Matrix prs = Matrix::allocate(nrows, (( terms.size()+pad-1 )/ pad ) * pad, 0); 
	//std::cout<<"in initialize_matrix"<<std::endl;
	//std::cout<<"seq_terms.size(): "<<seq_terms.size()<<" seq_coeffs.size(): "<<seq_coeffs.size()<<std::endl;
	//std::cout<<"terms.size(): "<<terms.size()<<std::endl;
	seq_coeffs_type::const_iterator cur_coeffs = seq_coeffs.begin();
   
	size_t i = 0;
   
	for(seq_terms_type::const_iterator cur_terms = seq_terms.begin();
			cur_coeffs != seq_coeffs.end();
			cur_coeffs++, cur_terms++, i++)
	{
		//std::cout<<"processing poly index: "<<i<<std::endl;
	   
	   assert(cur_terms->size() == cur_coeffs->size());
	   terms_type::const_iterator cur_term = (*cur_terms).begin();
	   //std::cout<<"cur_terms.size(): "<<cur_terms->size()<<" cur_coeffs.size(): "<<cur_coeffs->size()<<std::endl;
	   
	   for(coeffs_type::const_iterator cur_coeff = (*cur_coeffs).begin(); 
	       (cur_coeff != cur_coeffs->end()) && (cur_term != cur_terms->end()); 
	       cur_coeff++, cur_term++  )
	     {
		//     std::cout << "cur_term: " << *cur_term << std::endl ; 
		//     std::cout << "cur_coeff: " << *cur_coeff << std::endl;
		
		   size_t j = 0;		
		   for(terms_type::const_iterator it = terms.begin(); (it != terms.end())  ; it++, j++)
		     {		
			     //std::cout << "j: " << j << ", term: " << *it << std::endl;
			if(*cur_term == *it)
			{
			   pmatrix.setEntry(i, j, *cur_coeff);
			   break;			
			}	
		     }
		   
	     }
	   
		   
	}
}

void PolynomialSequence::order_pivot_ops(
		std::vector<F4Operations> &ops, std::vector<size_t> &deps) const
{
	typedef std::map<Term, uint32_t, Term::comparator> ordered_pivots_type;
	//std::map<Term, uint32_t, Term::comparator> pivotsOrdered(term_comparator);

	ordered_pivots_type pivotsOrdered( pivots.begin(), pivots.end(),
			term_comparator );
	//ordered_pivots_type pivotsOrdered( pivots.begin(), pivots.end(),
	//		term_comparator );
	//pivotsOrdered.insert(pivots.begin(), pivots.end());
	std::cout<<"in order_pivot_ops, pivotOpsOrdered.size(): "<<pivotsOrdered.size()<<std::endl;
	ops.push_back( F4Operations() );
	deps.assign(seq.size(), 0);

	#if 1 
	std::vector<size_t> l(seq.size(), 0);
	for(ordered_pivots_type::reverse_iterator it = pivotsOrdered.rbegin();
			it != pivotsOrdered.rend(); it++)
	{
		uint32_t o = it->second;
		std::cout<<"o: "<<o<<" monom: "<<it->first<<std::endl;
		// pivot_ops might not have it->first as a key
		pivot_ops_type::const_iterator entry = pivot_ops.find(it->first);;
		if (entry == pivot_ops.end()) {
			continue;
		}

		for(pivot_op_type::const_iterator j = entry->second.begin();
				j != entry->second.end(); j++)
		{
			uint32_t t = j->first;
			std::cout<<"t: "<<t<<"l[o]: "<<l[o]<<" l[t]: "<<l[t]<<" j->sec: "<<j->second<<std::endl;
			if(l[ o ] > l[t]) {
				l[t] = l[o];
			}
			ops[ l[t] ].push_back( t, o, j->second);
			deps[t]++; // target t has deps[t] operations do be done before it can be used as operator
			l[t]++; // one operation per level, attention this also affects the following if-statements

			if(l[t] >= ops.size()) {
				ops.push_back( F4Operations() );
			}
		}
	}
	#else
	size_t l = 0;
	for(std::map<const Term*, vector<std::pair<size_t, coeffType> >, TermComparator>::reverse_iterator it = pivotOpsOrdered.rbegin(); it != pivotOpsOrdered.rend(); it++)
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
	ops.pop_back();
	//prepareTime += seconds() - timer;
	std::cout << "Operations:\t" << ops.size() << "\n";
	//std::cout << "Preparation:\t" << seconds() - timer << "\n";
}

void PolynomialSequence::gauss(Matrix& pmatrix, std::vector<bool>& empty) const
{
	pmatrix.gauss(2*n_spolys, empty, field);
}

void PolynomialSequence::pReduce(Matrix& pmatrix,
		std::vector<F4Operations>& ops) const
{
	typedef std::vector<F4Operations> ops_type;
	for(ops_type::const_iterator i = ops.begin(); i != ops.end(); i++)
	{
		//		#pragma omp parallel for num_threads( threads ) 
		for(size_t j = 0; j < i->size(); j++) {
			pmatrix.pReduce(i->targets[j], i->opers[j],
					i->factors[j], field);
		}
	}
}

