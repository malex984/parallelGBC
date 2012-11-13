#ifndef F4MATRIX_H
#define F4MATRIX_H

#include <set>
#include <vector>
#include <map>

#include "CoeffField.H"



typedef std::vector<coeffType> CRow;
typedef std::vector<CRow*> BaseMatrix; 


class Matrix
{
  private:
    typedef BaseMatrix Base; 
    typedef coeffType Element;

    BaseMatrix* rows;

  public:

    Matrix(const int r, const int c);
    ~Matrix();

    Element getEntry(const int i, const int j) const;

    void setEntry(const int i, const int j, const Element& v);

    std::size_t size() const;
    size_t size(size_t row) const;


    /// elementary op. on rows: [target] -= [oper] * (factor)
    void pReduce(const std::size_t target, const std::size_t oper, const coeffType factor, const CoeffField* const field);
    
    /// empty: flags indicating zero rows
    void gauss(std::size_t upper, std::vector<bool>& empty, const CoeffField* const field);
}; // class Matrix

std::ostream& operator<< (std::ostream& out, const Matrix &M);

#endif // F4MATRIX_H
