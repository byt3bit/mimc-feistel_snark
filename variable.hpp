#ifndef VARIABLE_HPP_
#define VARIABLE_HPP_
#include "common/common.hpp"
#include "common/utility.hpp"

namespace snark {

  typedef uint32_t index_t;

  // template< typename field_t >
  // class svar;

  // template< typename field_t >
  // class svar_array;


  template< typename field_t > 
  class linear_combination;

  
  class variable
  {

  public:
    index_t index;
    variable(const index_t index = 0):
      index(index)
    {
    }
    
  };

  template< typename field_t >  
  class linear_term : public variable
  {

  public:
    field_t coef;
    
    linear_term() {};
    linear_term(const long unsigned int index, const field_t coef):
      variable(index), coef(coef)
    {
    }
    
    linear_term(const linear_term< field_t > &other);
    
    void set(const long unsigned int index, const field_t coef) {
      
      this->index = index; 
      this->coef = coef;
    }

    /* this computes A*(B_i*x_i) */
    linear_term< field_t > operator*(const field_t field_coef);
   
    friend std::ostream& operator<< (std::ostream &out, const linear_term &lt) {

      if(std::is_same<field_t, NTL::GF2>::value) {
	
	if(lt.index != 0) out<<" X_"<<lt.index;
	else
	  out<<" 1";
      } 
      else if(std::is_same<field_t, NTL::GF2E>::value){
	
	if(lt.index == 0) out<<" const";
	else
	  out<<lt.coef<<" *X_"<<lt.index;
      }
      return out;
    }

  };

  template< typename field_t >
  inline bool operator==(const linear_term< field_t > &lhs, const linear_term< field_t > &rhs) {
    
      return lhs.index == rhs.index && lhs.coef == rhs.coef ;
  }

  /* For further use */
  /* variable linear_term , allows to assign 'value' to the linear term */
  
  // template< typename field_t>
  // class linear_term_t: public linear_term< field_t > 
  // {

  // public:
  //   field_t val;

  //   linear_term_t(){};
  //   linear_term_t(const long unsigned int index, const field_t coef, const field_t val):
  //     linear_term< field_t > (index, coef), val(val)
  //   {
  //   }
  //   void set(const long unsigned int index, const field_t coef, const field_t val) {

  //     this->index = index;
  //     this->coef = coef;
  //     this->val = val;
  //   }

  //   void set(const long unsigned int index, const field_t val) {

  //     this->index = index;
  //     this->val = val;
  //   }

  // };
  
 
  /* data structure to represent the rank-1 constraints 
     of the for (A,w)*(B,w) = (C,w)

   */

  template< typename field_t > 
  class linear_combination {

  public:
    
    std::vector< linear_term< field_t > > lin_terms;
    
    linear_combination();

    linear_combination(linear_term< field_t > lt);

    linear_combination(const linear_combination< field_t > &other);

    /*constructed by a vector of linear terms*/
    linear_combination(const std::vector< linear_term< field_t > > &lts);

    linear_combination(std::initializer_list< field_t > lts); 
    
    typename std::vector<linear_term< field_t > >::const_iterator begin() const;
    typename std::vector<linear_term< field_t > >::const_iterator end() const;

    void clear() {
      
      if(!lin_terms.empty()) lin_terms.clear();
    }

    linear_term< field_t >& operator[](const size_t i); 
    
    size_t size()
    {
      return lin_terms.size();
    }
    
    void add_term(const linear_term< field_t > &lterm); 

    void reset(linear_combination < field_t > &other); 


    linear_combination< field_t > operator+(const linear_combination< field_t > &other_lc);
    
    //linear_combination< field_t > operator+(const linear_term< field_t > &lt);
    
    bool operator==(linear_combination< field_t > &other);

    field_t evaluate(std::vector< field_t > &var_assignment);

    friend std::ostream& operator<< (std::ostream &out, const linear_combination< field_t > &lc)
    {
      out<<"( ";
      size_t i = 0;
      for(i = 0; i < lc.lin_terms.size()-1;i++)
	out<<lc.lin_terms[i]<<" + ";
      out<<lc.lin_terms[i]<<" )";
      return out;
    }
 
  };

  template<typename field_t > 
  linear_combination< field_t > operator+( const linear_combination< field_t > &lc1, const linear_term< field_t > lt2 );

  template<typename field_t> 
  linear_combination< field_t > operator+( const linear_term< field_t > lt1, const linear_term< field_t > lt2 );

#include "variable.cpp"

}//mimc_snark

#endif
