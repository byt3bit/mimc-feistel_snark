#ifndef MIMC_CRF_SNARK_HH
#define MIMC_CRF_SNARK_HH
#include "common/common.hpp"
#include "r1_constraint.hpp"
#include "variable.hpp"

namespace snark {

  void generate_erf_roundconst(int r, int nbr, std::vector< std::vector< NTL::GF2E > > &rc);

  template<typename field_t> 
  class mimcfstl_crf_snark {

  private:

    unsigned int addCount;
    unsigned int multCount;

    uint16_t numofBranch;
    int numofRound;
    uint16_t branchSize;
    std::vector< field_t > roundConst;
        
  public:
    
    snarkcs< field_t > mimcfstlcrf_constr_wit;

    mimcfstl_crf_snark();
    mimcfstl_crf_snark(uint16_t numof_branch, uint16_t branch_size,
		       std::vector< field_t > &round_const):
    numofBranch(numof_branch), branchSize(branch_size), roundConst(round_const)
    {
      
      addCount = 0;multCount = 0;
      numofRound = (int) (2*(branchSize)/log2(3.0)) + 6*numofBranch - 5;

      assert(numofRound<= (int)roundConst.size());

    }


    void generate_r1_constraint();
    void generate_witness(std::vector<field_t>  ptxt);

    unsigned int num_of_add() {
      return addCount;
    }
    int num_of_round() {
      return numofRound;
    }

    unsigned int num_of_mult() {
      return multCount;
    }
  };
 

#include "mimccrf_snark.cc"

}


#endif
