#ifndef MIMC_FSTL_2BR_SNARK_HH
#define MIMC_FSTL_2BR_SNARK_HH
#include "common/common.hpp"
#include "r1_constraint.hpp"
#include "variable.hpp"

namespace snark {

  template<typename field_t>
  class mimcfstl_2br_snark {

  private:
    unsigned int addCount;
    unsigned int multCount;

    int numofRound;
    uint16_t branchSize;
    std::vector<field_t> roundConst;

  public:

    snarkcs< field_t > mimcfstl2br_constr_wit;

    mimcfstl_2br_snark();

    mimcfstl_2br_snark(uint16_t branch_size,
		   std::vector<field_t> &round_const):
      addCount(0), multCount(0), branchSize(branch_size), roundConst(round_const)
    
    {
      numofRound = (int) ((2*branchSize)/log2(3.0));
    }

    void generate_r1_constraint();
    void generate_witness(field_t leftinput, field_t rightinput);

    unsigned int num_of_add() {
      return addCount;
    }

    unsigned int num_of_mult() {
      return multCount;
    }


  };

#include "mimcfstl_2br_snark.cpp"

}


#endif
