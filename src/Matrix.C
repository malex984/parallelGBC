#include "../include/Matrix.H"

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
