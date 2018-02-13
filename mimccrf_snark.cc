#ifndef MIMC_CRF_SNARK_CPP
#define MIMC_CRF_SNARK_CPP


template<typename field_t> 
mimcfstl_crf_snark<field_t>::mimcfstl_crf_snark():
  addCount(0), multCount(0)
{
}


/*
This function computes all the rank-1 constraints for the MiMC Fesitel with 
CRF having t branches.
*/


template<typename field_t> 
void mimcfstl_crf_snark< field_t >::generate_r1_constraint()
{

  std::vector< uint16_t > indxary(numofBranch);

  index_t var_index = 1;
  
  while(var_index <= numofBranch) {
    
    indxary[(int)var_index - 1] = var_index;
    var_index++;
  }

  for(int jround = 0;jround < numofRound;jround++) {

    linear_term< field_t > xc(0, roundConst[jround]); 
    linear_combination< field_t > A(xc);

    index_t tmp = indxary[numofBranch - 1];
    for(int i = numofBranch - 2;i >= 0;i--) {
      
      linear_term< field_t > x( indxary[i], (field_t) 1 );
      A.add_term(x);
      indxary[i + 1] = indxary[i];
    }
    
    linear_combination< field_t > B(A);
    linear_term< field_t >  xt(var_index, (field_t) 1 );
    var_index++;

    linear_combination < field_t > C(xt);
    constraint < field_t > constr( A, B, C);
    mimcfstlcrf_constr_wit.add_constraint(constr);

    indxary[0] = var_index;

    A.reset(C);
    C.clear();
    xt.set(var_index, (field_t) 1);
    var_index++;

    C.add_term(xt);
    xt.set(tmp, (field_t) (-1) );
    C.add_term(xt);

    
    constr.reset_constraint(A, B, C);

    mimcfstlcrf_constr_wit.add_constraint(constr);
    
  }
  //mimcfstlcrf_constr_wit.print_constraints();

}

/*

This function computes the values of all the intermediate snark variables  
which are used to define the set rank-1 of constraints.                    
+--------------------------------------------------------------------------+
*/
 

template<typename field_t> 
void mimcfstl_crf_snark<field_t>::generate_witness(std::vector<field_t> input)
{

  index_t var_index = 0;
  std::vector<field_t> roundInput(input);
  
  while(var_index < numofBranch) { 

    mimcfstlcrf_constr_wit.add_witness( roundInput.at(var_index) );
    var_index++;
  }

  for(int jround = 0;jround < numofRound;jround++) {

    field_t temp;
    
    temp = roundConst[jround];

    for(int i = 1;i < numofBranch;i++) {
      temp += roundInput.at(i); 
      addCount++;
    }
    field_t yval = temp*temp;
    //var_index++;
    multCount++;
    mimcfstlcrf_constr_wit.add_witness(yval);

    temp = temp*yval;
    multCount++;
  
    temp = temp + roundInput.at(0);
    addCount++;
    //var_index++;
    mimcfstlcrf_constr_wit.add_witness(temp); 
    
    for( int i = 0;i < (int)roundInput.size()-1;i++) 
      roundInput.at(i) = roundInput.at(i+1); 

    roundInput.at(roundInput.size()-1) = temp;
  }
}






#endif
