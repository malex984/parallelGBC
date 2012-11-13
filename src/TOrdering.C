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
#include "TOrdering.H"
#include "Term.H"

  int DegRevLexOrdering::cmp(const Term& a, const Term& b) const
  {   
	  //std::cout<<"in DegRevLex::cmp, a: "<<a<<" b: "<<b<<std::endl;
    if(a == b) return 0;
    //std::cout<<"not equal"<<std::endl;
    if(a.deg() == b.deg())
    {   
	    //std::cout<<"degrees equal"<<std::endl;
      for(size_t i = N; i > 0; i--)
      {   
	      //std::cout<<"i: "<<i<<std::endl;
        degreeType r = a[i-1] - b[i-1];
	//std::cout<<"r: "<<r<<std::endl;
        if(r != 0)
        {   
          return r > 0 ? -1 : 1 ; 
        }   
      }   
      return 0;
    }   
    return a.deg() < b.deg() ? -1 : 1;
  }

  int LexOrdering::cmp(const Term& a, const Term& b) const 
  {
	  //std::cout<<"in Lex::cmp, a: "<<a<<" b: "<<b<<std::endl;
	if(a == b) return 0;
    //std::cout<<"not equal"<<std::endl;
	for(long i = 0; i < N; i++)
	{
	      //std::cout<<"i: "<<i<<std::endl;
		degreeType r = a[i] - b[i];
	//std::cout<<"r: "<<r<<std::endl;
		if(r != 0) {
			return r < 0 ? -1 : 1;
		}
	}
	return 0;
  }
