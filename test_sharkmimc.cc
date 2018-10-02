#include "common/ntl_header.hpp"
#include "common/utility.hpp"
#include "mimc_snark.hh"
#include "sharkmimc.hh"
#include "sharkmimc_invrf.hh"
#include <gmp.h>
#include <chrono>


int sharkmimc_round(int nbranch, int branch_size)
{

  int round_fullhat = ((nbranch + 2 < 2*branch_size) ? 6 : 8);
  int r = (int)ceil( log2(nbranch*branch_size + branch_size - 2) / log2((branch_size + 1)*0.5) ) ;
  int r1 =(int)ceil( log2(nbranch*branch_size + branch_size - 2) );
  
  int round_temp =  (int) (2*ceil( (r + r1)/2 ))  + 4;
  //n > 4 for snark which is branchSize
  int round_full = fmax(round_fullhat, round_temp);
  int round_partial = fmax(fmax(0, 0.63*branch_size + 4 - round_full), 0.63*(branch_size/2 + log2(nbranch)) + 2 - round_full);
  int nround = round_full + round_partial;
  //assert(numofRound <= (int) roundConst.size());
  assert(round_partial > 0);

  std::cout<<"num of round with full Sbox "<<round_full<<"\n";
  std::cout<<"num of round with partial Sbox "<<round_partial<<"\n";
  std::cout<<"num of rounds = "<<nround<<"\n";
  return nround;

}

int sharkmimc_invrf_round(int nbranch, int branch_size)
{

  int round_full = (((nbranch + 2) < branch_size) ? 6 : 8);
      
  int r = (int) ceil( 2*branch_size/log2((float) nbranch) );
  r = r + 4;
  r = r - round_full;
  int round_partial = (int) fmax(0, r);
  
  int nround = round_full + round_partial;

  std::cout<<"num of round with full Sbox "<<round_full<<"\n";
  std::cout<<"num of round with partial Sbox "<<round_partial<<"\n";
  std::cout<<"num of rounds = "<<nround<<"\n";

  return nround;

}

void sharkmimc_profile()
{


  int branch_size = 181;
  NTL::ZZ zp = NTL::GenPrime_ZZ(branch_size, 80);
  NTL::ZZ_p::init(zp);

  std::vector< std::vector< NTL::ZZ_p> > mds_Mat;
  std::vector< NTL::ZZ_p > input;
  std::vector< NTL::ZZ_p > round_const;
  int nbranch = 6;
  mds_Mat = snark::generate_mds_gfp(nbranch);
  int nround = sharkmimc_round(nbranch, branch_size);
  snark::generate_mimc_roundconst_gfp(round_const, nround*nbranch);

  for(int i = 0; i < branch_size;i++)
    input.push_back( NTL::random_ZZ_p() );

  float time_wit = 0.0, time_total = 0.0;
  float time_constr = 0.0; 
  int niter = 1000;
  for(int iter = 0; iter < niter;iter++) {

    auto time_start = std::chrono::high_resolution_clock::now();
    snark::sharkmimc_snark< NTL::ZZ_p > g_sharkmimc_snark(nbranch, branch_size, round_const);

    g_sharkmimc_snark.generate_r1_constraint(mds_Mat);
    auto time_end1 = std::chrono::high_resolution_clock::now();

    g_sharkmimc_snark.generate_witness(mds_Mat, input);
    auto time_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_fp_ms_constr = time_end1-time_start;

    std::chrono::duration<double, std::milli> time_fp_ms_wit = time_end-time_end1;
    std::chrono::duration<double, std::milli> time_fp_ms_total = time_end-time_start;
    
    time_constr += time_fp_ms_constr.count();
    time_wit += time_fp_ms_wit.count();
    time_total += time_fp_ms_total.count();
  }
  time_constr /= (float)niter;
  time_wit /= (float) niter;
  time_total /= (float) niter;
  std::cout<<"time to generate constraint "<<time_constr<<" ms\n";
  std::cout<<"time to generate witness "<<time_wit<<" ms\n";
  std::cout<<"time taken in total "<<time_total<<" ms\n";

}

void sharkmimc_invrf_profile() {

  int branch_size = 171;
  NTL::ZZ zp = NTL::GenPrime_ZZ(branch_size, 80);
  NTL::ZZ_p::init(zp);

  std::vector< std::vector< NTL::ZZ_p> > mds_Mat;
  std::vector< NTL::ZZ_p > round_const;

  std::vector< NTL::ZZ_p > input;
  int nbranch = 6;
  mds_Mat = snark::generate_mds_gfp(nbranch);
  int nround = sharkmimc_invrf_round(nbranch, branch_size);
  snark::generate_mimc_roundconst_gfp(round_const, nround*nbranch);

  for(int i = 0; i < branch_size;i++)
    input.push_back( NTL::random_ZZ_p() );

  float time_wit = 0.0, time_total = 0.0;
  float time_constr = 0.0; 
  int niter = 2000;
  for(int iter = 0; iter < niter;iter++) {

    auto time_start = std::chrono::high_resolution_clock::now();
    snark::sharkmimc_invrf_snark< NTL::ZZ_p > g_sharkmimc_invrf_snark(nbranch, branch_size, round_const);

    g_sharkmimc_invrf_snark.generate_r1_constraint(mds_Mat);
    auto time_end1 = std::chrono::high_resolution_clock::now();
    
    g_sharkmimc_invrf_snark.generate_witness(mds_Mat, input);
    auto time_end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double, std::milli> time_fp_ms_constr = time_end1-time_start;

    std::chrono::duration<double, std::milli> time_fp_ms_wit = time_end-time_end1;
    std::chrono::duration<double, std::milli> time_fp_ms_total = time_end-time_start;
    
    time_constr += time_fp_ms_constr.count();
    time_wit += time_fp_ms_wit.count();
    time_total += time_fp_ms_total.count();
  }
  time_constr /= (float)niter;
  time_wit /= (float) niter;
  time_total /= (float) niter;
  
  std::cout<<"time to generate constraint "<<time_constr<<" ms\n";
  std::cout<<"time to generate witness "<<time_wit<<" ms\n";
  std::cout<<"time taken in total "<<time_total<<" ms\n";

}


int main() {


  //sharkmimc_profile();
  sharkmimc_invrf_profile();

  //NTL::ZZ_p rc = NTL::random_ZZ_p();
  

  NTL::vec_ZZ_p a;
  a.SetLength(5);

  // mpz_t gprim, rp;
  // std::stringstream ssa;
  // ssa<<zp;
  // mpz_set_str(gprim, ssa.str().c_str(), 10);
  // mpz_cdiv_q_ui (rp, gprim, (unsigned long)2);
  // //gmp_printf("ap = %Zd \n", ap);
  // char bignumstr[50];
  // mpz_get_str (bignumstr, 10, rp);
  // NTL::ZZ_p u;
  // NTL::conv(u, bignumstr);
  // //std::cout<<"bp = "<<u<<"\n";
  // std::cout<<"rc = "<<rc<<"\n";
  // rc = (NTL::ZZ_p) 1;
  // std::cout<<"rc = "<<rc<<"\n";

  
  return 0;
}
