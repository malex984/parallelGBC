#include "Matrix.H"

using namespace std;

Matrix::Matrix(const int r, const int c)
{
	rows = new BaseMatrix();
	std::cout<<"r: "<<r<<" c: "<<c<<std::endl;
	//resize(r);
	for(int i = 0; i < r; i++) {
		std::cout<<"i: "<<i<<std::endl;
		rows->push_back(new CRow(c, 0));
	}
	std::cout<<"size: "<<rows->size()<<std::endl;
	//assign(r, CRow(c, 0));
}

Matrix::~Matrix()
{
	for(size_t i = 0; i < rows->size(); i++) {
		delete rows->at(i);
	}
}

Matrix::Element Matrix::getEntry(const int i, const int j) const
{
	return rows->at(i)->at(j);
}

void Matrix::setEntry(const int i, const int j, const Element& v)
{
rows->at(i)->at(j) = v;
}

std::size_t Matrix::size() const
{
return rows->size();
}

size_t Matrix::size(size_t row) const
{
return rows->at(row)->size();
}

void Matrix::pReduce(const std::size_t target, const std::size_t oper, const coeffType factor, const CoeffField* const field)
{
	std::cerr << "pReduce(target: " << target << ", oper: " << oper << ", factor: " << factor << std::endl;

	const coeffType modn = field->getChar();
	const coeffType c = field->getLog( factor );

	std::cerr << "c: " << c << std::endl;
	std::cerr << "modn: " << modn << std::endl;

	for(std::size_t k = 0; k < this->size(oper); k++)
	{
		std::cerr << "[" << k  << "]: " ;

		const coeffType ok = this->getEntry(oper, k);
		std::cerr << "ok: " << ok;

		if(ok != 0)
		{
			//          const coeffType b = field->getExp(field->getLog(ok) + c);

			const coeffType tk = this->getEntry(target, k);
			//          std::cerr << ", tk: " << tk;
			//          std::cerr << "b: " << b << std::endl;
			//          this->setEntry(target, k, (b > tk) ? tk - b + modn : tk - b);

			//          std::cerr << "tk': " << this->getEntry(target, k) << std::endl;
			//          std::cerr << "tk'': " << field->mulSub(tk, ok, c) << std::endl;
			this->setEntry(target, k, field->mulSub(tk, ok, c));

			std::cerr << "tk': " << this->getEntry(target, k) << ", must be: " << field->mulSub(tk, ok, c);
		}

		std::cerr << std::endl;
	}
}

    /// empty: flags indicating zero rows
void Matrix::gauss(std::size_t upper, std::vector<bool>& empty, const CoeffField* const field)
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


ostream& operator<< (ostream& out, const Matrix &M)
{
   out << "[" << endl;
   
   for(size_t i = 0; i < M.size(); i++)
   {
      out << "\t/" << i << "/:\t"<< M.getEntry(i, 0);
      
      for(size_t j = 1; j < M.size(i); j++)
	out << ",\t"<< M.getEntry(i, j);
      
      out << endl;
   }
   
   out << "]";
   return out;
}
