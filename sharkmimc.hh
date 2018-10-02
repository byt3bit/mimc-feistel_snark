#ifndef SHARKMIMC_SNARK_HH
#define SHARKMIMC_SNARK_HH


namespace snark {

  std::vector< std::vector < NTL::ZZ_p> > generate_mds_gfp(int t);

  template<typename field_t> 
  class sharkmimc_snark {

    private:

    unsigned int addCount;
    unsigned int multCount;

    uint16_t numofBranch;
    int numofRound;
    uint16_t branchSize;
    std::vector< field_t >  roundConst;

    int round_full;
    int round_partial;

    public:
    
    snarkcs< field_t > sharkmimc_constr_wit;

    sharkmimc_snark();
    sharkmimc_snark(uint16_t numof_branch, uint16_t branch_size, std::vector < field_t > & round_const):
      addCount(0), multCount(0), numofBranch(numof_branch), branchSize(branch_size), roundConst(round_const)
    {
      int round_fullhat = ((numofBranch + 2 < 2*branchSize) ? 6 : 8);
      int round_temp = (int)ceil(log2(numofBranch*branchSize + branchSize -2)) + (int)(ceil(log2(numofBranch*branchSize + branchSize -2)/log2((branchSize + 1)/2))) + 4;
      //n > 4 for snark which is branchSize
      round_full = fmax(round_fullhat, round_temp);
      round_partial = fmax(fmax(0, 0.63*branchSize + 4 - round_full), 0.63*(branchSize/2 + log2(numofBranch)) + 2 - round_full);
      numofRound = round_full + round_partial;
      //assert(numofRound <= (int) roundConst.size());
      assert(round_partial > 0);

      // std::cout<<"num of round with full Sbox "<<round_full<<"\n";
      // std::cout<<"num of round with partial Sbox "<<round_partial<<"\n";
      //std::cout<<"num of rounds = "<<numofRound<<"\n";

    }

    void generate_r1_constraint( std::vector< std::vector< field_t> >  mds_Mat );
    void generate_witness( std::vector< std::vector< field_t> >  mds_Mat, std::vector<field_t> ptxt);

    unsigned int num_of_add() {
      return addCount;
    }
    // int num_of_round() {
    //   return numofRound;
    // }

    unsigned int num_of_mult() {
      return multCount;
    }

    int num_of_round() {
      return numofRound;
    }


  };
  
#include "sharkmimc.cc"
}


#endif
