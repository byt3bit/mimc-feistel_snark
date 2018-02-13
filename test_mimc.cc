#include "variable.hpp"
#include "common/ntl_header.hpp"
#include "mimc_snark.hh"
#include <chrono>
#include <iomanip>


void profile_mimc_snark_gf2e(int blockSize, int niter, bool ver= true);

void profile_mimc_snark_gfp(int blockSize, int niter, bool ver=true);


int main() 
{


  int num_iter = 2000;
  //average time over num_iter iterations
  
  //profile_mimc_snark_gf2e(1024, numiter);
  //over GF(2^n)
  
  profile_mimc_snark_gfp(1024, num_iter, true);
  //over GF(p)

}


void profile_mimc_snark_gf2e(int blockSize, int niter, bool ver) {

  NTL::GF2X modp = NTL::BuildIrred_GF2X(blockSize);
  NTL::GF2E::init(modp);
  std::vector< NTL::GF2E > round_const;
  
  int nround = (int) (blockSize/log2(3.0));

  snark::generate_mimc_roundconst(round_const, blockSize, nround);
    
  assert(round_const.size() == (unsigned int)nround);

  float time_wit = 0.0, time_constr = 0.0, time_total = 0.0;
  int nmult(0), nadd(0);
  

  for(int i = 0; i < niter; i++) {
    
    NTL::GF2E z = NTL::random_GF2E();  /*random input*/

    auto time_start = std::chrono::high_resolution_clock::now();

    snark::mimc_em_snark<NTL::GF2E> g_mimcp_snark(blockSize, round_const);
    
  
    g_mimcp_snark.generate_r1_constraint();
  
    auto time_end1 = std::chrono::high_resolution_clock::now();

    g_mimcp_snark.generate_witness(z);

    auto time_end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> time_fp_ms_constr = time_end1-time_start;
    std::chrono::duration<double, std::milli> time_fp_ms_wit = time_end-time_end1;
    std::chrono::duration<double, std::milli> time_fp_ms_total = time_end-time_start;

    time_wit += time_fp_ms_wit.count();
    time_constr += time_fp_ms_constr.count();
    time_total += time_fp_ms_total.count();

    if(!i) {
      nmult = g_mimcp_snark.num_of_mult();
      nadd = g_mimcp_snark.num_of_addition();
    }
    
  }

  time_wit /= (float)niter;
  time_constr /= (float)niter;
  time_total /= (float)niter;

  if(ver) {

    printf("\nMiMC profile: \n\n");

    std::cout<<std::left<<std::setw(20)<<"#rounds "<<nround<<"\n";
    std::cout<<std::setw(20)<<"#Multiplications "<<nmult<<"\n";
    std::cout<<std::setw(20)<<"#Additions "<<nadd<<"\n";
    std::cout<<std::setw(35)<<FXD_MSG1<<std::setprecision(5)<<time_constr<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG2<<std::setprecision(5)<<time_wit<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG3<<std::setprecision(5)<<time_total<<" ms\n";
    

    
  }
}

void profile_mimc_snark_gfp(int blockSize, int niter, bool ver) {

  NTL::ZZ zp = NTL::GenPrime_ZZ(blockSize, 80);
  NTL::ZZ_p::init(zp);
  std::vector< NTL::ZZ_p > round_const;
  
  int nround = (int) (blockSize/log2(3.0));

  
  snark::generate_mimc_roundconst_gfp(round_const, nround);
  
  assert(round_const.size() == (size_t)nround);
  int nmult(0), nadd(0);

  float time_wit = 0.0, time_constr = 0.0, time_total = 0.0;

  for(int i = 0; i < niter; i++) {
    
    NTL::ZZ_p z = NTL::random_ZZ_p();  /*random input*/

    auto time_start = std::chrono::high_resolution_clock::now();
    
    snark::mimc_em_snark< NTL::ZZ_p > g_mimcp_snark(blockSize, round_const);
  
    g_mimcp_snark.generate_r1_constraint();
  
    auto time_end1 = std::chrono::high_resolution_clock::now();

    g_mimcp_snark.generate_witness(z);

    auto time_end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> time_fp_ms_constr = time_end1-time_start;
    std::chrono::duration<double, std::milli> time_fp_ms_wit = time_end-time_end1;
    std::chrono::duration<double, std::milli> time_fp_ms_total = time_end-time_start;

    time_wit += time_fp_ms_wit.count();
    time_constr += time_fp_ms_constr.count();
    time_total += time_fp_ms_total.count();

    if(!i) {
      nmult = g_mimcp_snark.num_of_mult();
      nadd = g_mimcp_snark.num_of_addition();
    }
    
  }

  time_wit /= (float)niter;
  time_constr /= (float)niter;
  time_total /= (float)niter;

  if(ver) {

    printf("\nMiMC profile (over F_p): \n\n");
    std::cout<<std::left<<std::setw(20)<<"Block size "<<blockSize<<"\n";
    
    std::cout<<std::setw(20)<<"#rounds "<<nround<<"\n";
    std::cout<<std::setw(20)<<"#Multiplications "<<nmult<<"\n";
    std::cout<<std::setw(20)<<"#Additions "<<nadd<<"\n";
    std::cout<<std::setw(35)<<FXD_MSG1<<std::setprecision(5)<<time_constr<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG2<<std::setprecision(5)<<time_wit<<" ms\n";
    std::cout<<std::setw(35)<<FXD_MSG3<<std::setprecision(5)<<time_total<<" ms\n";

    
  }
}



