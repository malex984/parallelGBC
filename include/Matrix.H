#ifndef F4MATRIX_H
#define F4MATRIX_H

#include <set>
#include <vector>
#include <map>

#include "../include/CoeffField.H"


#ifdef HAVE_LELA

# include "lela/util/commentator.h"
#include "lela/blas/context.h"
// #include "lela/ring/gf2.h"
#include "lela/ring/modular.h"

#include "lela/blas/level3.h"
#include "lela/solutions/echelon-form.h"
#include "lela/solutions/echelon-form-gf2.h"

#endif



#ifdef HAVE_LELA

typedef LELA::Modular<coeffType> LELARing;
typedef LELA::DenseMatrix<typename LELARing::Element> LelaDenseMatrix;
typedef LELA::SparseMatrix<typename LELARing::Element> LelaSparseMatrix;


using namespace LELA;

/*
EchelonForm<Ring>::METHOD_STANDARD_GJ;
EchelonForm<Ring>::METHOD_ASYMPTOTICALLY_FAST_GJ;
EchelonForm<Ring>::METHOD_FAUGERE_LACHARTRE;
EchelonForm<Ring>::METHOD_UNKNOWN;

EchelonForm<GF2>::METHOD_M4RI
*/

template <class Ring, class LELAMatrix>
void my_row_echelon_form (const Ring &R, LELAMatrix& A,
                          typename EchelonForm<Ring>::Method method,
                          bool reduced)
{
  Context<Ring> ctx (R);

  commentator.start ("Converting matrix to row-echelon-form", __FUNCTION__);

  if (method == EchelonForm<Ring>::METHOD_FAUGERE_LACHARTRE && !reduced) {
    reduced = true;
    commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
        << "Note: constructing reduced row-echelon form, since non-reduced form is not available with the chosen method" << std::endl;
  }

  EchelonForm<Ring> EF (ctx);

  if (reduced)
    commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
        << "Computing reduced form" << std::endl;
  else
    commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
        << "Computing non-reduced form" << std::endl;

  try
  {
    EF.echelonize (A, reduced, method);
  }
  catch (LELAError e)
  {
    commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR) << e;
    commentator.stop ("error");
    throw e;
  }
  catch (...)
  {
    commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
        << "Non-LELA exception occured during 'echelonize'" << std::endl;
    commentator.stop ("error");
    throw;
  }

  commentator.stop (MSG_DONE);
}

#endif





#ifdef HAVE_LELA
typedef LelaDenseMatrix BaseMatrix; 
// typedef LelaSparseMatrix BaseMatrix;
#else
typedef std::vector<coeffType> CRow;
typedef vector<CRow> BaseMatrix; 
#endif


class Matrix: private BaseMatrix  
{
  private:
    typedef BaseMatrix Base; 
#ifndef HAVE_LELA
    typedef coeffType Element;
#else
    typedef typename Base::Element Element;
#endif

  public:
    Matrix(): Base ()
    {
    }

#ifdef HAVE_LELA
    Matrix(const int r, const int c): Base(r, c)
#else
    Matrix(const int r, const int c) // : Base(r, CRow(c))
#endif
    {
      assign(r, CRow(c, 0));
    }

   // TODO: what about destructor???

    Element getEntry(const int i, const int j) const
    {
#ifndef HAVE_LELA
      return (*(Base*)this)[i][j];
#else
      Element x = 0;
      return ((Base*)this)->getEntry (x, i, j);
#endif
    }

    void setEntry(const int i, const int j, const Element& v)
    {
#ifndef HAVE_LELA
      (*(Base*)this)[i][j] = v;
#else
//      typename LelaDenseMatrix::Element x; R.init(x, v);
      cerr << "setEntry[" << i << ", " << j << ", [" << v << "]) => ";
      ((Base*)this)->setEntry(i, j, v);
      cerr << getEntry(i, j) << endl;
#endif
    }

    std::size_t size() const
    {
#ifndef HAVE_LELA
      return ((Base*)this)->size();
#else
      return rowdim();
#endif
    }

    std::size_t size(std::size_t row) const
    {
#ifndef HAVE_LELA
      return (*(Base*)this)[row].size();
#else
      return coldim();
#endif
    }

    static Matrix* allocate(const int r, const int c, const coeffType& defv = 0)
    {
#ifndef HAVE_LELA
      Matrix* p = new Matrix(); p->assign(r, CRow(c, defv) ); return p;
#else
      return new Matrix(r, c);      
#endif
    }

    /// elementary op. on rows: [target] -= [oper] * (factor)
    void pReduce(const std::size_t target, const std::size_t oper, const coeffType factor, const CoeffField* const field)
    {
//      cerr << "pReduce(target: " << target << ", oper: " << oper << ", factor: " << factor << endl;
       
//      const coeffType modn = field->getChar();
      const coeffType c = field->getLog( factor );

//      cerr << "c: " << c << endl;
//      cerr << "modn: " << modn << endl;
      
      for(std::size_t k = 0; k < this->size(oper); k++)
      {
//        cerr << "[" << k  << "]: " ;
	 
        const coeffType ok = this->getEntry(oper, k);
//        cerr << "ok: " << ok;
	 
        if(ok != 0)
        {
//          const coeffType b = field->getExp(field->getLog(ok) + c);

	  const coeffType tk = this->getEntry(target, k);
//          cerr << ", tk: " << tk;
//          cerr << "b: " << b << endl;
//          this->setEntry(target, k, (b > tk) ? tk - b + modn : tk - b);
	   
//          cerr << "tk': " << this->getEntry(target, k) << endl;
//          cerr << "tk'': " << field->mulSub(tk, ok, c) << endl;
          this->setEntry(target, k, field->mulSub(tk, ok, c));

//          cerr << "tk': " << this->getEntry(target, k) << ", must be: " << field->mulSub(tk, ok, c);
        }
	 
        cerr << endl;
      }
    }

    /// empty: flags indicating zero rows
    void gauss(std::size_t upper, std::vector<bool>& empty, const CoeffField* const field)
    {
//      LELARing R(field->getChar());
//      LELA::my_row_echelon_form(R, *this, LELA::EchelonForm<LELARing>::METHOD_UNKNOWN, true);
      
      for(std::size_t i = 1; i < upper; i+=2)
      {
        std::size_t p = 0;
        bool found = false;
        coeffType factor = 0;
        
        for(p = 0; !found && p < size(i); p++) {
          factor = getEntry(i, p);
          found  = (factor != 0);
        }
        p--;
        empty[i] = !found;
        
        if(found)
        {
      // Normalize
          if(factor != 1)
          {
            factor = field->inv(factor);
            for(size_t j = p; j < size(i); j++)
            {
              setEntry(i, j, field->mul(getEntry(i, j), factor));
            }
          }
      // Execute
//			#pragma omp parallel for num_threads( threads )
          for(std::size_t j = 2; j < upper; j+=2)
          {
            const std::size_t k = (i+j)%upper;
            
            if( getEntry(k, p) != 0) {
              coeffType factor = field->getLog(getEntry(k, p));
              for(std::size_t m = p; m < size(k); m++)
              {
            // This is mulSub for primitives not vectors !
                setEntry(k, m, field->mulSub(getEntry(k, m), getEntry(i, m), factor));
              }
            }
          }
        }
      }
    }
}; // class Matrix

ostream& operator<< (ostream& out, const Matrix &M);

#endif // F4MATRIX_H
