#ifndef MIMC_SNARKBOARD_HPP_
#define MIMC_SNARKBOARD_HPP_

#include "variable.hpp"
#include "common/common.hpp"

namespace snark {


  template<typename field_t>
  class constraint {
 
  public: 
    linear_combination< field_t > A, B, C;
    
    constraint();
    constraint(const linear_combination<field_t> a,
	       const linear_combination<field_t> b,
	       const linear_combination<field_t> c):
      A(a), B(b), C(c)
    { 
    }
    constraint(const constraint &other):
      A(other.A), B(other.B), C(other.C)
    {
      //set_constraint_type(other.get_constraint_type());
    }

    void reset_constraint(const linear_combination<field_t> a,
			  const linear_combination<field_t> b,
			  const linear_combination<field_t> c);


    void clear();

    friend std::ostream& operator<< (std::ostream &out, const constraint &constr)
    {
      out<<constr.A<<" * "<<constr.B<<" = "<<constr.C;
      return out;
    }

  };
  
  template<typename field_t>
  class constraint_system {

  public:
    std::vector<constraint<field_t> > constraints;

    constraint_system() {}; 
    constraint_system(const constraint_system &other):
      constraints(other.constraints)
    {
    } 

    void add_constraint(const constraint<field_t> constr);
    constraint<field_t>& operator[](const size_t i);    
    
  };



  template<typename field_t> 
  class snarkcs {
    
  public:

    constraint_system<field_t> r1_constraints;
    std::vector<field_t> witness;

    snarkcs();
    snarkcs(const snarkcs &other):
      r1_constraints(other.r1_constraints), witness(other.witness)
    {
    }

    const size_t constraint_size() const;
    const size_t witness_size() const;

    void add_constraint(const constraint<field_t> constr); 
    void add_witness(const field_t w);
    void add_witness(std::initializer_list< field_t > l);
    void print_constraints() const; 
    void print_witness() const;

    /*removes the last constraint in the system of constraints*/
    void pop_constraint();
  };

#include "r1_constraint.cpp"

}

#endif
