#ifndef SHARKMIMC_INVRF_SNARK_CC
#define SHARKMIMC_INVRF_SNARK_CC

template<typename field_t> 
sharkmimc_invrf_snark<field_t>::sharkmimc_invrf_snark():
  addCount(0), multCount(0)
{
}

template<typename field_t> 
void sharkmimc_invrf_snark< field_t >::generate_r1_constraint( std::vector< std::vector< field_t> >  mds_Mat ) {

  std::vector< long unsigned int > indxary(numofBranch);

  index_t var_index = 1;

  while(var_index <= numofBranch) {
    
    indxary[var_index - 1] = var_index;
    var_index++;
    
  }

  int roundfull_first = (int)ceil(round_full/2);
  //int roundfull_last = round_full - roundfull_first;
  
  int jround = 0;
  std::vector<long unsigned int> temp_indx_arr(numofBranch);
  for(int jround = 0;jround < roundfull_first;jround++) {
    
    if(jround == 0) {
      
      for(int ibranch = 0;ibranch < numofBranch;ibranch++) {
      
	linear_term < field_t > xc(0, roundConst[jround*numofBranch + ibranch]);
	
	linear_term < field_t > x(indxary[numofBranch-1], (field_t) 1 );
	linear_combination< field_t > A({xc, x});
	linear_term < field_t > xt(var_index, (field_t) 1);
	linear_combination < field_t > B(xt);
	temp_indx_arr[ibranch] = var_index;
	var_index++;
	xc.set(0, (field_t) 1);
	linear_combination < field_t > C(xc);
    
	constraint < field_t > constr(A, B, C);
	sharkmimc_constr_wit.add_constraint(constr); // X * Y = 1
		
      }
    }
    else {
            
      memset(&indxary, 0, indxary.size() * sizeof(indxary[0]));
      indxary = temp_indx_arr;
      
      for(int i = 0; i < numofBranch;i++) {

	linear_combination <field_t> A;
	linear_term < field_t > xc(0, roundConst[jround*numofBranch + i]);
	A.add_term(xc);
	for(int jcol = 0;jcol < numofBranch;jcol++) {
	  linear_term<field_t> x (indxary[numofBranch - jcol - 1], mds_Mat[i][jcol]);
	  A.add_term(x);
	}
	linear_term <field_t> y(var_index, (field_t) 1);
	linear_combination <field_t> B(y);
	temp_indx_arr[i] = var_index;
	var_index++;
	xc.set(0, (field_t) 1);
	linear_combination <field_t> C(xc);
	constraint < field_t > constr(A, B, C);
	sharkmimc_constr_wit.add_constraint(constr); // X * Y = 1

      } 
    }
    
  } //jround ends

  //* rounds in the middle with partial Sbox  
  int ipartial = 0;
  for(;ipartial < round_partial;ipartial++) {
    
    memset(&indxary, 0, indxary.size() * sizeof(indxary[0]));
    indxary = temp_indx_arr;
   
    linear_combination <field_t> A;
    linear_term < field_t > xc(0, roundConst[jround*numofBranch + 0]);
    A.add_term(xc);
    for(int jcol = 0;jcol < numofBranch;jcol++) {
      linear_term<field_t> x (indxary.at(numofBranch - jcol - 1), mds_Mat[0][jcol]);
      A.add_term(x);
    }
    linear_term <field_t> y(var_index, (field_t) 1);
    linear_combination <field_t> B(y);
    temp_indx_arr[0] = var_index;
    var_index++;
    xc.set(0, (field_t) 1);
    linear_combination <field_t> C(xc);
    constraint < field_t > constr(A, B, C);
    sharkmimc_constr_wit.add_constraint(constr); // X * Y = 1
    
    //** This part is storing only the linear constraints and can be optimized (perhaps)
    
    for(int ibranch = 1;ibranch < numofBranch;ibranch++) {
      
      A.clear();
      xc.set(0, roundConst[jround*numofBranch + ibranch]);
      A.add_term(xc);
      for(int jcol = 0;jcol < numofBranch;jcol++) {
	linear_term<field_t> x (indxary.at(numofBranch - jcol - 1), mds_Mat[ibranch][jcol]);
	A.add_term(x);
      }
      xc.set(0, (field_t) 1);
      B.add_term(xc);
      linear_term <field_t> y(var_index, (field_t) 1);
      temp_indx_arr[ibranch] = var_index;
      var_index++;
      linear_combination <field_t> C(y);
      constraint < field_t > constr(A, B, C);
      sharkmimc_constr_wit.add_constraint(constr);
    }
    jround++;
    
  }//ipartial loop ends
  
  //* last full sbox rounds
  for(;jround < numofRound;jround++) {
    
    memset(&indxary, 0, indxary.size() * sizeof(indxary[0]));
    indxary = temp_indx_arr;
    
    for(int i = 0; i < numofBranch;i++) {
      
      linear_combination <field_t> A;
      linear_term < field_t > xc(0, roundConst[jround*numofBranch + i]);
      A.add_term(xc);
      for(int jcol = 0;jcol < numofBranch;jcol++) {
	linear_term<field_t> x (temp_indx_arr.at(numofBranch - jcol - 1), mds_Mat[i][jcol]);
	A.add_term(x);
      }
      linear_term <field_t> y(var_index, (field_t) 1);
      linear_combination <field_t> B(y);
      var_index++;
      xc.set(0, (field_t) 1);
      linear_combination <field_t> C(xc);
      constraint < field_t > constr(A, B, C);
      sharkmimc_constr_wit.add_constraint(constr); // X * Y = 1
      
    }
  }

}
template<typename field_t>
void sharkmimc_invrf_snark< field_t >::generate_witness( std::vector< std::vector< field_t> >  mds_Mat ,
						   std::vector< field_t > ptext)
{

  std::vector< field_t > temp(ptext);

  index_t var_index = 0;
  while(var_index < numofBranch) {

    sharkmimc_constr_wit.add_witness( temp.at(var_index) );
    var_index++;
  }

  int roundfull_first = (int)ceil(round_full/2);

  int jround = 0;
  
  for(;jround <  roundfull_first;jround++) {

    for(int ibranch = 0;ibranch < numofBranch;ibranch++) {
      
      temp[ibranch] = 1/(temp[ibranch] + roundConst[jround*numofBranch + ibranch]);
      sharkmimc_constr_wit.add_witness(temp[ibranch]);
      
    }
    std::vector< field_t > tvec(numofBranch);
    
    for(int i = 0; i < numofBranch;i++) {
      
      field_t r = (field_t) 0;
      for(int jcol = 0;jcol < numofBranch;jcol++) 
	r += temp[jcol]*mds_Mat[i][jcol];

      tvec[i] = r;
    }
    temp = tvec;
  }

  // * PARTIAL SBOX LAYER
  int ipartial = 0;
  for(;ipartial < round_partial;ipartial++) {

    temp[0] = 1/(temp[0] + roundConst[jround*numofBranch + 0]);

    sharkmimc_constr_wit.add_witness(temp[0]);

    for(int i = 1; i < numofBranch;i++) {
      
      temp[i] += roundConst[jround*numofBranch + i];
      sharkmimc_constr_wit.add_witness(temp[i]);
    }
    
    std::vector< field_t > tvec(numofBranch);
    for(int i = 0; i < numofBranch;i++) {
      
      field_t r = (field_t) 0;
      for(int jcol = 0;jcol < numofBranch;jcol++) 
	r += temp[jcol]*mds_Mat[i][jcol];

      tvec[i] = r;
    }
    temp = tvec;
    jround++;
  }//ipartial loop ends

  //* LAST FULL SBOX LAYER
  for(;jround < numofRound;jround++) {

    for(int ibranch = 0;ibranch < numofBranch;ibranch++) {
      
      temp[ibranch] = 1/(temp[ibranch] + roundConst[jround*numofBranch + ibranch]);
      sharkmimc_constr_wit.add_witness(temp[ibranch]);
      
    }
    std::vector< field_t > tvec(numofBranch);
    if(jround != numofRound - 1) {
      
      for(int i = 0; i < numofBranch;i++) {
      
	field_t r = (field_t) 0;
	for(int jcol = 0;jcol < numofBranch;jcol++) 
	  r += temp[jcol]*mds_Mat[i][jcol];

	tvec[i] = r;
      }
      temp = tvec;
    }
  }//jround loop ends

}





#endif
