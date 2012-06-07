
#ifndef POLYNOMIAL_SEQUENCE_H
#define POLYNOMIAL_SEQUENCE_H

#include "../include/Term.H"
#include "../include/Polynomial.H"
#include "../include/CoeffField.H"
#include "../include/Matrix.H"

#include <vector>
#include <utility>
#include <boost/unordered_map.hpp>

class F4Pivot {
	public:
		const Term* column;
		size_t row;
		vector<pair<size_t, coeffType> > operations;

		F4Pivot(const Term* column, size_t row) : column(column), row(row) {}

		void addOperation(size_t row, coeffType factor) {
			operations.push_back(make_pair(row, factor));
		}
};

class F4Operation {
	public:
		size_t target;
		size_t oper;
		coeffType factor;

		F4Operation(size_t target, size_t oper, coeffType factor) : target(target), oper(oper), factor(factor) {}
};


class PolynomialSequence {
public:
	typedef std::vector<std::pair<const Term*, Polynomial> > seq_type;
	typedef std::vector<coeffType> coeffs_type;
	typedef std::vector<coeffs_type > seq_coeffs_type;
	typedef std::vector<const Term*> terms_type;
	typedef std::vector<terms_type> seq_terms_type;;
	typedef boost::unordered_map<const Term*, size_t> pivots_type;
	typedef boost::unordered_map<const Term*, vector<pair<size_t, coeffType> > >
		pivot_ops_type;

	PolynomialSequence(degreeType deg, TMonoid& monoid,
			TermComparator& tog, CoeffField* field) :
		sugar_degree(deg), monoid(monoid), term_comparator(tog),
		field(field),
		seq(), seq_terms(), seq_coeffs(),terms(), pivots(),
		pivot_ops(), n_spolys(0) {}

	~PolynomialSequence() {}

	void push_spoly(const Term* lcm, const Polynomial& left,
			const Polynomial& right);
	void push_back(const Polynomial& poly);
	void push_back(const Term* t, const Polynomial& poly,
			bool check=true, bool pivot=true);
	template<typename ContainerType>
	void set_terms(typename ContainerType::const_iterator beg,
			typename ContainerType::const_iterator end) {
		terms.insert(terms.begin(), beg, end);
	};
	const coeffs_type& get_coeff_vector(size_t ind) const;
	const terms_type& get_terms_vector(size_t ind) const;
	void reduce(std::vector<Polynomial> &polys) const;

	const size_t size() const { return seq.size(); }
	bool check_pivot(const Term* t) const {
		return pivots.find(t) != pivots.end();
	}

protected:
	void initialize_matrix(Matrix& pmatrix) const;
	void get_new_polys(std::vector<Polynomial> &polys, Matrix& pmatrix,
			std::vector<bool> empty_rows) const;
	void order_pivot_ops(vector<vector<F4Operation> > &ops) const;
	void pReduce(Matrix& pmatrix, vector<vector<F4Operation> >& ops) const;
	void gauss(Matrix& pmatrix, vector<bool>& empty) const;
	

private:
	// this is passed to the constructor of polynomials
	degreeType sugar_degree;
	// monoid is used to construct new monomials
	const TMonoid& monoid;
	const TermComparator& term_comparator;
	CoeffField* field;
	seq_type seq;
	seq_terms_type seq_terms;
	seq_coeffs_type seq_coeffs;
	std::vector<const Term*> terms;
	pivots_type pivots;
	pivot_ops_type pivot_ops;
	size_t n_spolys;
};

#endif
