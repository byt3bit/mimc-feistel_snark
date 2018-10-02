#ifndef SHARKMIMC_SNARK_CC
#define SHARKMIMC_SNARK_CC


//** generate MDS (cauchy matrix) 
std::vector< std::vector < NTL::ZZ_p> >  generate_mds_gfp(int t) {

  // NTL::mat_ZZ_p mds_mat;
  // mds_mat.SetDims(t, t);
  NTL::vec_ZZ_p a, b;
  a.SetLength(t);
  b.SetLength(t);
  
  assert(NTL::IsZero(a));
  assert(NTL::IsZero(b));

  a[0] = NTL::random_ZZ_p();
  std::vector< std::vector < NTL::ZZ_p > > mdsmat;

  // if(!NTL::IsZero(a)) NTL::clear(a);
  // if(!NTL::IsZero(b)) NTL::clear(b);
  
  for(int i = 0;i < t;i++) {

    if(i) {
      a[i] = a[i-1] + (NTL::ZZ_p) 1;
    }
    if(!NTL::IsZero(b)) NTL::clear(b);
    std::vector< NTL::ZZ_p > temp_vec;
    for(int j = 0;j < t;j++) {
      if(j)
	b[j] = b[j-1] + 1;

      //mds_mat[i][j] = NTL::inv(a[i] + b[j]);
      temp_vec.push_back(NTL::inv(a[i] + b[j]));

    }
    mdsmat.push_back(temp_vec);
  }
  //NTL::mat_ZZ_p temp_mat = mds_mat;
  //assert(NTL::gauss(mds_mat) == t); //*check full rank/invertible.
  return mdsmat;
  
}

template<typename field_t> 
sharkmimc_snark<field_t>::sharkmimc_snark():
  addCount(0), multCount(0)
{
}

template<typename field_t> 
void sharkmimc_snark< field_t >::generate_r1_constraint( std::vector< std::vector< field_t> >  mds_Mat ) {

  std::vector< long unsigned int > indxary(numofBranch);

  index_t var_index = 1;

  while(var_index <= numofBranch) {
    
    indxary[var_index - 1] = var_index;
    var_index++;
    
  }

  int roundfull_first = (int)ceil(round_full/2);
  
  int jround = 0;
  std::vector<long unsigned int> temp_indx_arr(numofBranch);
  
  for(int jround = 0;jround < roundfull_first;jround++) {
    
    if(jround == 0) {
      
      for(int ibranch = 0;ibranch < numofBranch;ibranch++) {
      
	linear_term < field_t > xc(0, roundConst[jround*numofBranch + ibranch]);
	
	linear_term < field_t > x(indxary[ibranch], (field_t) 1 );
	linear_combination< field_t > A({xc, x});

	linear_combination < field_t > B(A);
	linear_term < field_t > xt(var_index, (field_t) 1);
	var_index++;
	linear_combination < field_t > C(xt);
    
	constraint < field_t > constr(A, B, C);
	sharkmimc_constr_wit.add_constraint(constr); // X * X = Y
	A.reset(C);
	xt.set(var_index, (field_t) 1);
	temp_indx_arr[ibranch] = var_index;
	var_index++;
	C.clear();
	C.add_term(xt);
	constr.reset_constraint(A, B, C);
	sharkmimc_constr_wit.add_constraint(constr); // X * Y = Z
		
      }
    }
    else {
      
      memset(&indxary, 0, indxary.size() * sizeof(indxary[0]));
      //**collect the variables which will be combined using mds matrix
      indxary = temp_indx_arr;
      
      for(int i = 0; i < numofBranch;i++) {

	linear_combination <field_t> A;
	linear_term < field_t > xc(0, roundConst[jround*numofBranch + i]);
	A.add_term(xc);
	for(int jcol = 0;jcol < numofBranch;jcol++) {
	  linear_term<field_t> x (indxary.at(numofBranch - jcol - 1), mds_Mat[i][jcol]);
	  A.add_term(x);
	}
	linear_combination <field_t> B(A);
	linear_term <field_t> y(var_index, (field_t) 1);
	var_index++;
	linear_combination <field_t> C(y);
	constraint < field_t > constr(A, B, C);
	sharkmimc_constr_wit.add_constraint(constr); // X * X = Y
	A.reset(C);
	y.set(var_index, (field_t) 1);
	temp_indx_arr[i] = var_index;
	var_index++;
	C.clear();
	C.add_term(y);
	constr.reset_constraint(A, B, C);
	sharkmimc_constr_wit.add_constraint(constr); // X * Y = Z
      } 
    }
    
  }

  //* rounds with partial Sbox in the middle 
  int ipartial = 0;
  for(;ipartial < round_partial;ipartial++) {
    
    memset(&indxary, 0, indxary.size() * sizeof(indxary[0]));
    //**collect the variables which will be combined using mds matrix
    indxary = temp_indx_arr;

    linear_combination <field_t> A;
    linear_term < field_t > xc(0, roundConst[jround*numofBranch + 0]);
    A.add_term(xc);
    for(int jcol = 0;jcol < numofBranch;jcol++) {
      linear_term<field_t> x (indxary.at(numofBranch - jcol - 1), mds_Mat[0][jcol]);
      A.add_term(x);
    }
    linear_combination <field_t> B(A);
    linear_term <field_t> y(var_index, (field_t) 1);
    var_index++;
    linear_combination <field_t> C(y);
    constraint < field_t > constr(A, B, C);
    sharkmimc_constr_wit.add_constraint(constr); // X * X = Y
    A.reset(C);
    y.set(var_index, (field_t) 1);
    temp_indx_arr[0] = var_index;
    var_index++;
    C.clear();
    C.add_term(y);
    constr.reset_constraint(A, B, C);
    sharkmimc_constr_wit.add_constraint(constr); // X * Y = Z

    //** This part is storing only the linear constraints and can be optimized (?)

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

  }

  //* LAST FULL SBOX ROUNDS
  for(;jround < numofRound;jround++) {

    memset(&indxary, 0, indxary.size() * sizeof(indxary[0]));
    //**collect the variables which will be combined using mds matrix
    indxary = temp_indx_arr;

    for(int i = 0; i < numofBranch;i++) {
      
      linear_combination <field_t> A;
      linear_term < field_t > xc(0, roundConst[jround*numofBranch + i]);
      A.add_term(xc);
      for(int jcol = 0;jcol < numofBranch;jcol++) {
	
	linear_term<field_t> x (indxary.at(numofBranch - jcol - 1), mds_Mat[i][jcol]);
	A.add_term(x);
      }
      linear_combination <field_t> B(A);
      linear_term <field_t> y(var_index, (field_t) 1);
      var_index++;
      linear_combination <field_t> C(y);
      constraint < field_t > constr(A, B, C);
      sharkmimc_constr_wit.add_constraint(constr); // X * X = Y
      A.reset(C);
      y.set(var_index, (field_t) 1);
      temp_indx_arr[i] = var_index;
      var_index++;
      C.clear();
      C.add_term(y);
      constr.reset_constraint(A, B, C);
      sharkmimc_constr_wit.add_constraint(constr); // X * Y = Z
    }
  }

}

template<typename field_t> 
void sharkmimc_snark< field_t >::generate_witness( std::vector< std::vector<field_t> >  mds_Mat ,
						   std::vector< field_t > ptext)
{

  std::vector< field_t > temp(ptext);
  int jround = 0;
  int roundfull_first = (int)ceil(round_full/2);

  index_t var_index = 0;

  while(var_index < numofBranch) {

    sharkmimc_constr_wit.add_witness( temp.at(var_index) );
    var_index++;
  }
  
  
  for(;jround < roundfull_first;jround++) {

    for(int ibranch = 0;ibranch < numofBranch;ibranch++) {
      
      field_t r = temp[ibranch] + roundConst[jround*numofBranch + ibranch];
      temp[ibranch] = r*r;
      multCount++;
      sharkmimc_constr_wit.add_witness(temp[ibranch]);
      temp[ibranch] *= r;
      multCount++;
      sharkmimc_constr_wit.add_witness(temp[ibranch]);
      
    }
    std::vector< field_t > tempv(numofBranch);
    for(int irow = 0; irow < numofBranch;irow++) {
      field_t r = (field_t) 0;
      for(int jcol = 0;jcol < numofBranch;jcol++) {

	r += mds_Mat[irow][jcol]*temp[jcol];
	if(jcol) addCount++;
      }
      tempv[irow] = r;
      
    }
    temp = tempv;
    
  }//jround loop ends

  int ipartial = 0;
  for(;ipartial < round_partial;ipartial++) {

    for(int ibranch = 0;ibranch < numofBranch;ibranch++)
      temp[ibranch] += roundConst[jround*numofBranch + ibranch];

    field_t r = temp[0];
    temp[0] = r*r;
    multCount++;
    sharkmimc_constr_wit.add_witness(temp[0]);
    temp[0] *= r;
    multCount++;
    sharkmimc_constr_wit.add_witness(temp[0]);

    for(int i = 1;i < numofBranch;i++)
      sharkmimc_constr_wit.add_witness(temp[i]);

    std::vector< field_t > tempv(numofBranch);
    for(int irow = 0; irow < numofBranch;irow++) {
      field_t r = (field_t) 0;
      for(int jcol = 0;jcol < numofBranch;jcol++) {

	r += mds_Mat[irow][jcol]*temp[jcol];
	if(jcol) addCount++;
      }
      tempv[irow] = r;
    }
    temp = tempv;

    jround++;
  }

  // ** LAST FULL SBOX LAYER
  for(;jround < numofRound;jround++) {

    for(int ibranch = 0;ibranch < numofBranch;ibranch++) {
      
      field_t r = temp[ibranch] + roundConst[jround*numofBranch + ibranch];
      
      temp[ibranch] = r*r;
      multCount++;
      sharkmimc_constr_wit.add_witness(temp[ibranch]);
      temp[ibranch] *= r;
      multCount++;
      sharkmimc_constr_wit.add_witness(temp[ibranch]);
    }
    
    std::vector< field_t > tempv(numofBranch);
    if(jround != numofRound - 1 ) {
      
      for(int irow = 0; irow < numofBranch;irow++) {
	
	field_t r = (field_t) 0;
	for(int jcol = 0;jcol < numofBranch;jcol++) {

	  r += mds_Mat[irow][jcol]*temp[jcol];
	  if(jcol) addCount++;
	}
	tempv[irow] = r;
      }
      temp = tempv;
    }

  }

}





#endif
