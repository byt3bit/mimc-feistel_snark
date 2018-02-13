#ifndef MIMC_ERF_SNARK_CPP
#define MIMC_ERF_SNARK_CPP

void generate_erf_roundconst(int r, int nbr, int br_sz, std::vector< std::vector< NTL::GF2E > > &rc)
{
  //generate constant when the structure is defined with distinct round round constants for each
  //output of the ER function

  for(int i = 0;i < r;i++) {

    std::vector< NTL::GF2E > tvec;
    for(int t = 0; t < nbr - 1 ;t++) { 
      
      NTL::GF2X tmp;
      NTL::GF2E c;

      for(int j = 0;j < br_sz;j++) {
	NTL::SetCoeff(tmp, (long)j, (NTL::GF2)getrandbit());
      }
      NTL::conv(c, tmp); //conv to GF2E from GF2X
      tvec.emplace_back(c);
    }
    rc.emplace_back(tvec);
    
  }
  
}

template<typename field_t> 
mimcfstl_erf_snark<field_t>::mimcfstl_erf_snark():
  addCount(0), multCount(0)
{
}


template<typename field_t> 
void mimcfstl_erf_snark< field_t >::generate_r1_constraint() {

  std::vector< long unsigned int > indxary(numofBranch);

  index_t var_index = 1;

  while(var_index <= numofBranch) {
    
    indxary[var_index - 1] = var_index;
    var_index++;
  }

  for(int jround = 0;jround < numofRound;jround++) {

    linear_term < field_t > xc(0, roundConst[jround]); 
    linear_term < field_t > x(indxary[numofBranch-1], (field_t) 1 );
    linear_combination< field_t > A({xc, x});

    linear_combination < field_t > B(A);
    linear_term < field_t > xt(var_index, (field_t) 1);
    var_index++;
    linear_combination < field_t > C(xt);
    

    constraint < field_t > constr(A, B, C);
    mimcfstlerf_constr_wit.add_constraint(constr);
    A.reset(C);

    long unsigned int tempi = indxary[numofBranch - 1];
    for(int i = numofBranch - 2;i >= 0;i--) {

      x.set(var_index, (field_t) 1 );
      C.clear();
      C.add_term(x);
      indxary[i + 1] = var_index;
      var_index++;
      
      x.set(indxary[i], (field_t)(-1) );
      C.add_term(x);
      
      constr.reset_constraint(A, B, C);
      mimcfstlerf_constr_wit.add_constraint(constr);

    }
    indxary[0] = tempi;
     
  } //jorund loop ends
  
  //mimcfstlerf_constr_wit.print_constraints();
 
}

/*
// Generates the witness vector for the MiMC-erf 
*/

template<typename field_t> 
void mimcfstl_erf_snark<field_t>::generate_witness(std::vector<field_t>& input) {

  index_t var_index = 0;
  std::vector<field_t> roundInput(input);
  
  while(var_index < numofBranch) {
    
    mimcfstlerf_constr_wit.add_witness( roundInput.at(var_index) );
    var_index++;
  }

  for(int jround = 0;jround < numofRound;jround++) {


    field_t temp = roundConst[jround];
    
    temp += roundInput[numofBranch - 1];
    field_t tmp1 = temp*temp;
    mimcfstlerf_constr_wit.add_witness( tmp1 );
    temp = temp*tmp1;
    mimcfstlerf_constr_wit.add_witness( temp );
    multCount += 2;
      
    field_t tmp2 = roundInput[numofBranch - 1];
    for(int i = numofBranch - 2;i >= 0;i--) {
      
      tmp1 = temp + roundInput.at(i);
      addCount++;
      roundInput[i + 1] = tmp1;
      mimcfstlerf_constr_wit.add_witness( tmp1 ); //causing more time for erf ??
      // numofBranch-1 witness storing causing more time for erf ??
    }
    roundInput[0] = tmp2;
    
  }//jround loop ends

}




#endif
