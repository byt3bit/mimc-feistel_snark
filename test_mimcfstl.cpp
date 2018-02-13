#include <vector>
#include <stdint.h>
#include <chrono>
#include <ctime>
#include <iomanip>
#include "common/ntl_header.hpp"
#include "common/utility.hpp"
#include "mimc_snark.hh"
#include "mimcerf_snark.hh"
#include "mimccrf_snark.hh"
#include "mimcfstl_2br_snark.hh"


void profile_mimcfstl_crf_gf2e(uint16_t numofBranch, uint16_t branchSize, 
			  int numofIter, bool ver=true) {


  NTL::GF2X modp = NTL::BuildIrred_GF2X(branchSize);
  NTL::GF2E::init(modp);

  std::vector< NTL::GF2E > ptext;
  std::vector< NTL::GF2E > roundConst;
  
  int nround = (int) (2*branchSize/log2(3.0)) + 6*numofBranch - 5;

  snark::generate_mimc_roundconst(roundConst, branchSize, nround);

  for(int i = 0;i < numofBranch;i++) {     /* random input */
    
    NTL::GF2E z = NTL::random_GF2E();
    ptext.emplace_back(z);
    
  }

  assert(roundConst.size() == (size_t) nround);
  
  float time_wit = 0.0, time_constr = 0.0, time_total = 0.0;

  unsigned int nmult, nadd;
  
  for(int iter = 0;iter < numofIter;iter++) {

    auto time_start = std::chrono::high_resolution_clock::now();


    snark::mimcfstl_crf_snark<NTL::GF2E> g_mimcfstl_snark(numofBranch, branchSize, roundConst);

    g_mimcfstl_snark.generate_r1_constraint();

    auto time_end1 = std::chrono::high_resolution_clock::now();
    
    g_mimcfstl_snark.generate_witness(ptext);

    auto time_end = std::chrono::high_resolution_clock::now();
    
  
    std::chrono::duration<double, std::milli> time_fp_ms_constr = time_end1-time_start;
    std::chrono::duration<double, std::milli> time_fp_ms_wit = time_end-time_end1;
    std::chrono::duration<double, std::milli> time_ms_total = time_end-time_start;


    time_constr += time_fp_ms_constr.count();
    time_wit += time_fp_ms_wit.count();
    time_total += time_ms_total.count();

    if(!iter) {
      
      nadd = g_mimcfstl_snark.num_of_add();
      nmult = g_mimcfstl_snark.num_of_mult();
    }
    
  }//iter loop ends

  
  time_wit /= (float)numofIter;
  time_constr /= (float)numofIter;
  time_total /=(float)numofIter;

  if(ver) {

    printf("\nMiMC Feistel-CRF profile: \n\n");

    std::cout<<std::left<<std::setw(30)<<"No. of branches "<<numofBranch<<"\n";
    std::cout<<std::setw(30)<<"Branch size "<<branchSize<<"\n\n";

    std::cout<<std::left<<std::setw(35)<<"#multiplications "<<nmult<<"\n";
    std::cout<<std::setw(35)<<"#additions "<<nadd<<"\n";
    std::cout<<std::setw(35)<<"#rounds "<<nround<<"\n";
    std::cout<<std::setw(35)<<FXD_MSG1<<std::setprecision(5)<<time_constr<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG2<<std::setprecision(5)<<time_wit<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG3<<std::setprecision(5)<<time_total<<" ms\n";

  }

}


void profile_mimcfstl_crf_gfp(uint16_t numofBranch, uint16_t branchSize, 
			      int numofIter, bool ver=true) {



  NTL::ZZ zp = NTL::GenPrime_ZZ(branchSize, 80);
  NTL::ZZ_p::init(zp);
  std::vector< NTL::ZZ_p > round_const;

  std::vector< NTL::ZZ_p > ptext;
  std::vector< NTL::ZZ_p > roundConst;
  
  int nround = (int) (2*branchSize/log2(3.0)) + 6*numofBranch - 5;

  snark::generate_mimc_roundconst_gfp(roundConst, nround);

  assert(roundConst.size() == (size_t) nround);


  float time_wit = 0.0, time_constr = 0.0, time_total = 0.0;

  unsigned int nmult, nadd;
  
  for(int iter = 0;iter < numofIter;iter++) {

    for(int i = 0;i < numofBranch;i++) {     /* random input */
    
      NTL::ZZ_p z = NTL::random_ZZ_p();
      ptext.emplace_back(z);
    
    }

    auto time_start = std::chrono::high_resolution_clock::now();


    snark::mimcfstl_crf_snark< NTL::ZZ_p > g_mimcfstl_snark(numofBranch, branchSize, roundConst);

    g_mimcfstl_snark.generate_r1_constraint();

    auto time_end1 = std::chrono::high_resolution_clock::now();
    
    g_mimcfstl_snark.generate_witness(ptext);

    auto time_end = std::chrono::high_resolution_clock::now();
    
  
    std::chrono::duration<double, std::milli> time_fp_ms_constr = time_end1-time_start;
    std::chrono::duration<double, std::milli> time_fp_ms_wit = time_end-time_end1;
    std::chrono::duration<double, std::milli> time_ms_total = time_end-time_start;


    time_constr += time_fp_ms_constr.count();
    time_wit += time_fp_ms_wit.count();
    time_total += time_ms_total.count();

    if(!iter) {
      
      nadd = g_mimcfstl_snark.num_of_add();
      nmult = g_mimcfstl_snark.num_of_mult();
      //nround = g_mimcfstl_snark.num_of_round();
    }
    
  }//iter loop ends

  
  time_wit /= (float)numofIter;
  time_constr /= (float)numofIter;
  time_total /=(float)numofIter;

  if(ver) {

    printf("\nMiMC Feistel-CRF profile (over F_p): \n\n");

    std::cout<<std::left<<std::setw(30)<<"No. of branches "<<numofBranch<<"\n";
    std::cout<<std::setw(30)<<"Branch size "<<branchSize<<"\n\n";

    std::cout<<std::setw(35)<<"#multiplications "<<nmult<<"\n";
    std::cout<<std::setw(35)<<"#additions "<<nadd<<"\n";
    std::cout<<std::setw(35)<<"#rounds "<<nround<<"\n";
    std::cout<<std::setw(35)<<FXD_MSG1<<std::setprecision(5)<<time_constr<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG2<<std::setprecision(5)<<time_wit<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG3<<std::setprecision(5)<<time_total<<" ms\n";

  }

}


void profile_mimcfstl_2br(uint16_t branchSize, int numofIter, bool ver=true) {


  int numofBranch = 2;
  
  NTL::GF2X modp = NTL::BuildIrred_GF2X(branchSize);
  NTL::GF2E::init(modp);

  std::vector< NTL::GF2E > roundConst;
  int nround = (int) (2*branchSize/log2(3.0));

  snark::generate_mimc_roundconst(roundConst, branchSize, nround);

  NTL::GF2E leftinput, rightinput;

  /* random input */
  leftinput = NTL::random_GF2E();
  rightinput = NTL::random_GF2E();
  
  assert(roundConst.size() == (size_t) nround);

  std::cout<<std::left<<std::setw(30)<<"No. of branches "<<numofBranch<<"\n";
  std::cout<<std::setw(30)<<"Branch size "<<branchSize<<"\n";
  std::cout<<std::setw(30)<<"Number of rounds "<<nround<<"\n\n";


  float time_wit = 0.0, time_constr = 0.0, time_total = 0.0;

  unsigned int nmult, nadd;
  
  for(int iter = 0;iter < numofIter;iter++) {

    auto time_start = std::chrono::high_resolution_clock::now();
    
    snark::mimcfstl_2br_snark<NTL::GF2E> g_mimcfstl_snark(branchSize, roundConst);

    g_mimcfstl_snark.generate_r1_constraint();

    auto time_end1 = std::chrono::high_resolution_clock::now();
    
    g_mimcfstl_snark.generate_witness(leftinput, rightinput);

    auto time_end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> time_fp_ms_constr = time_end1-time_start;
    std::chrono::duration<double, std::milli> time_fp_ms_wit = time_end-time_end1;
    std::chrono::duration<double, std::milli> time_ms_total = time_end-time_start;

    time_constr += time_fp_ms_constr.count();
    time_wit += time_fp_ms_wit.count();
    time_total += time_ms_total.count();
    
    if(!iter) {
      nadd = g_mimcfstl_snark.num_of_add();
      nmult = g_mimcfstl_snark.num_of_mult();
    }

    

  } //iter loop ends
  time_constr /= numofIter;
  time_wit /= numofIter;
  time_total /= numofIter;
  

  if(ver) {

    printf("MiMC 2 branch Feistel profile: \n\n");

    std::cout<<std::left<<std::setw(35)<<"#multiplications "<<nmult<<"\n";
    std::cout<<std::setw(35)<<"#additions "<<nadd<<"\n";
    std::cout<<std::setw(35)<<"#rounds "<<nround<<"\n";
    std::cout<<std::setw(35)<<FXD_MSG1<<std::setprecision(5)<<time_constr<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG2<<std::setprecision(5)<<time_wit<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG3<<std::setprecision(5)<<time_total<<" ms\n";
  }
  
}


void profile_mimcfstl_2br_gfp(uint16_t branchSize, int numofIter, bool ver=true) {

  NTL::ZZ zp = NTL::GenPrime_ZZ(branchSize, 80);
  NTL::ZZ_p::init(zp);
  std::vector< NTL::ZZ_p > roundConst;
  
  int nround = (int) (2*branchSize/log2(3.0));

  snark::generate_mimc_roundconst_gfp(roundConst, nround);

  NTL::ZZ_p leftinput, rightinput;

  /* random input */
  leftinput = NTL::random_ZZ_p();
  rightinput = NTL::random_ZZ_p();
  
  assert(roundConst.size() == (size_t) nround);

  float time_wit = 0.0, time_constr = 0.0, time_total = 0.0;

  unsigned int nmult, nadd;
  
  for(int iter = 0;iter < numofIter;iter++) {

    auto time_start = std::chrono::high_resolution_clock::now();
    
    snark::mimcfstl_2br_snark< NTL::ZZ_p > g_mimcfstl_snark(branchSize, roundConst);

    g_mimcfstl_snark.generate_r1_constraint();

    auto time_end1 = std::chrono::high_resolution_clock::now();
    
    g_mimcfstl_snark.generate_witness(leftinput, rightinput);

    auto time_end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> time_fp_ms_constr = time_end1-time_start;
    std::chrono::duration<double, std::milli> time_fp_ms_wit = time_end-time_end1;
    std::chrono::duration<double, std::milli> time_ms_total = time_end-time_start;

    time_constr += time_fp_ms_constr.count();
    time_wit += time_fp_ms_wit.count();
    time_total += time_ms_total.count();
    
    if(!iter) {
      nadd = g_mimcfstl_snark.num_of_add();
      nmult = g_mimcfstl_snark.num_of_mult();
      //nround = g_mimcfstl_snark.num_of_round();
    }

    

  } //iter loop ends
  time_constr /= numofIter;
  time_wit /= numofIter;
  time_total /= numofIter;
  

  if(ver) {

    printf("MiMC 2 branch Feistel profile (over F_p): \n\n");

    
    std::cout<<std::left<<std::setw(30)<<"Branch size "<<branchSize<<"\n";
    
    std::cout<<std::setw(35)<<"#multiplications "<<nmult<<"\n";
    std::cout<<std::setw(35)<<"#additions "<<nadd<<"\n";
    std::cout<<std::setw(35)<<"#rounds "<<nround<<"\n";
    std::cout<<std::setw(35)<<FXD_MSG1<<std::setprecision(5)<<time_constr<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG2<<std::setprecision(5)<<time_wit<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG3<<std::setprecision(5)<<time_total<<" ms\n";
  }
  
}

void profile_mimcfstl_erf(uint16_t numofBranch, uint16_t branchSize, 
			  int numofIter, bool ver=true) {

  NTL::GF2X modp = NTL::BuildIrred_GF2X(branchSize);
  NTL::GF2E::init(modp);

  std::vector< NTL::GF2E > ptext;
  std::vector< NTL::GF2E > roundConst;
  int nround  =  (int) (2*(branchSize)/log2(3.0)) + 4*numofBranch - 3;

  snark::generate_mimc_roundconst(roundConst, branchSize, nround);

  
  for(int i = 0;i < numofBranch;i++) {     /* random input */
    
    NTL::GF2E z = NTL::random_GF2E();
    ptext.emplace_back(z);
    
  }

  
  assert(roundConst.size() == (size_t) nround);

   
  float time_wit = 0.0, time_constr = 0.0, time_total = 0.0;

  unsigned int nmult, nadd;
  
  for(int iter = 0;iter < numofIter;iter++) {

    auto time_start = std::chrono::high_resolution_clock::now();


    snark::mimcfstl_erf_snark<NTL::GF2E> g_mimcfstl_snark(numofBranch, branchSize, roundConst);

    g_mimcfstl_snark.generate_r1_constraint();

    auto time_end1 = std::chrono::high_resolution_clock::now();
    
    g_mimcfstl_snark.generate_witness(ptext);

    auto time_end = std::chrono::high_resolution_clock::now();
    
  
    std::chrono::duration<double, std::milli> time_fp_ms_constr = time_end1-time_start;
    std::chrono::duration<double, std::milli> time_fp_ms_wit = time_end-time_end1;
    std::chrono::duration<double, std::milli> time_ms_total = time_end-time_start;


    time_constr += time_fp_ms_constr.count();
    time_wit += time_fp_ms_wit.count();
    time_total += time_ms_total.count();
    
    if(!iter) {
      
      nadd = g_mimcfstl_snark.num_of_add();
      nmult = g_mimcfstl_snark.num_of_mult();
    }
    
  }//iter loop ends

  
  time_wit /= (float)numofIter;
  time_constr /= (float)numofIter;
  time_total /= (float)numofIter;

  if(ver) {

    printf("\nMiMC Feistel-ERF profile: \n\n");

    std::cout<<std::left<<std::setw(30)<<"No. of branches "<<numofBranch<<"\n";
    std::cout<<std::setw(30)<<"Branch size "<<branchSize<<"\n\n";


    std::cout<<std::left<<std::setw(35)<<"#multiplications "<<nmult<<"\n";
    std::cout<<std::setw(35)<<"#additions "<<nadd<<"\n";
    std::cout<<std::setw(35)<<"#rounds "<<nround<<"\n";
    std::cout<<std::setw(35)<<FXD_MSG1<<std::setprecision(5)<<time_constr<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG2<<std::setprecision(5)<<time_wit<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG3<<std::setprecision(5)<<time_total<<" ms\n";

  }


}



void profile_mimcfstl_erf_gfp(uint16_t numofBranch, uint16_t branchSize, 
			  int numofIter, bool ver=true) {


  NTL::ZZ zp = NTL::GenPrime_ZZ(branchSize, 80);
  NTL::ZZ_p::init(zp);
  
  std::vector< NTL::ZZ_p > ptext;
  std::vector< NTL::ZZ_p > roundConst;
  int nround  =  (int) (2*(branchSize)/log2(3.0)) + 4*numofBranch - 3;

  snark::generate_mimc_roundconst_gfp(roundConst, nround);
  
  for(int i = 0;i < numofBranch;i++) {     /* random input */
    
    NTL::ZZ_p z = NTL::random_ZZ_p();
    ptext.emplace_back(z);
    
  }

  
  assert(roundConst.size() == (size_t) nround);

   
  float time_wit = 0.0, time_constr = 0.0, time_total = 0.0;

  unsigned int nmult, nadd;
  
  for(int iter = 0;iter < numofIter;iter++) {

    auto time_start = std::chrono::high_resolution_clock::now();


    snark::mimcfstl_erf_snark<NTL::ZZ_p> g_mimcfstl_snark(numofBranch, branchSize, roundConst);

    g_mimcfstl_snark.generate_r1_constraint();

    auto time_end1 = std::chrono::high_resolution_clock::now();
    
    g_mimcfstl_snark.generate_witness(ptext);

    auto time_end = std::chrono::high_resolution_clock::now();
    
  
    std::chrono::duration<double, std::milli> time_fp_ms_constr = time_end1-time_start;
    std::chrono::duration<double, std::milli> time_fp_ms_wit = time_end-time_end1;
    std::chrono::duration<double, std::milli> time_ms_total = time_end-time_start;


    time_constr += time_fp_ms_constr.count();
    time_wit += time_fp_ms_wit.count();
    time_total += time_ms_total.count();
    
    if(!iter) {
      
      nadd = g_mimcfstl_snark.num_of_add();
      nmult = g_mimcfstl_snark.num_of_mult();
    }
    
  }//iter loop ends

  
  time_wit /= (float)numofIter;
  time_constr /= (float)numofIter;
  time_total /= (float)numofIter;

  if(ver) {

    printf("\nMiMC Feistel-ERF profile (over F_p): \n\n");

    std::cout<<std::left<<std::setw(30)<<"No. of branches "<<numofBranch<<"\n";
    std::cout<<std::setw(30)<<"Branch size "<<branchSize<<"\n\n";


    std::cout<<std::left<<std::setw(35)<<"#multiplications "<<nmult<<"\n";
    std::cout<<std::setw(35)<<"#additions "<<nadd<<"\n";
    std::cout<<std::setw(35)<<"#rounds "<<nround<<"\n";
    std::cout<<std::setw(35)<<FXD_MSG1<<std::setprecision(5)<<time_constr<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG2<<std::setprecision(5)<<time_wit<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG3<<std::setprecision(5)<<time_total<<" ms\n";

  }


}





int main() {



  //profile_mimcfstl_crf_gf2e( 16, 65, 2000 );
  //printf("\n");
  //over GF(2^n)


  //profile_mimcfstl_crf_gfp(8, 128, 2000);
  //printf("\n");

  profile_mimcfstl_2br_gfp(512, 2000);
  printf("\n");


  profile_mimcfstl_erf_gfp(4, 256, 2000);
  printf("\n");
  profile_mimcfstl_erf_gfp(8, 128, 2000);
  printf("\n");
  profile_mimcfstl_erf_gfp(13, 80, 2000);
  printf("\n");

  // profile_mimcfstl_erf_gfp(16, 64, 2000);
  // printf("\n");*/

  return 0;
}
