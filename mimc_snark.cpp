#ifndef MIMC_SNARK_CPP
#define MIMC_SNARK_CPP

/*
//  generates round constants over F_{2^n}
*/

void generate_mimc_roundconst(std::vector< NTL::GF2E > &roundconst, 
			      int blocksize, int numround)
{

  for(int i = 0;i < numround;i++) {
    NTL::GF2X tmp;
    NTL::GF2E c;

    for(int j = 0;j < blocksize;j++) {
      NTL::SetCoeff(tmp, (long)j, (NTL::GF2)getrandbit());
    }
    NTL::conv(c, tmp); //conv to GF2E from GF2X
    roundconst.emplace_back(c);
  }
}

/*
//  generates round constants over F_p
*/

void generate_mimc_roundconst_gfp(std::vector< NTL::ZZ_p > &roundconst, 
			      int numround)
{

  for(int i = 0;i < numround;i++) {

    NTL::ZZ_p rc = NTL::random_ZZ_p();   
    roundconst.emplace_back(rc);

  }
  
}

template<typename field_t>
mimc_em_snark<field_t>::mimc_em_snark():
  xorCount(0), multCount(0) 
{
}

template<typename field_t> 
void mimc_em_snark<field_t>::generate_r1_constraint()
{

  index_t var_index = 1;
  
  for(int jround = 0;jround < num_round;jround++) {

    
    linear_term< field_t > x0(0, roundConst[jround]);
    linear_term < field_t > x1(var_index, (field_t) 1 );
    linear_term < field_t > y(var_index+1, (field_t) 1 );
    linear_term < field_t > z(var_index+2, (field_t) 1 );
    
    linear_combination< field_t > A(x0 + x1);
    linear_combination< field_t > B(A), C(y);


    //std::cout<<x1<<"\n";
    
    constraint<field_t> constr(A, B, C);  
    mimc_constr_wit.add_constraint(constr);

    A.reset(C);
    C.clear();
    C.add_term(z);
    constr.reset_constraint(A, B, C);
    mimc_constr_wit.add_constraint(constr);
    
    var_index += 2;
  }

  //mimc_constr_wit.print_constraints();

}


template<typename field_t> 
void mimc_em_snark<field_t>::generate_witness(field_t xval)
{

  //index_t var_index = 1;
  mimc_constr_wit.add_witness(xval);

  for(int jround = 0;jround < num_round;jround++) {

    field_t temp = xval + roundConst[jround];
    xorCount++;

    field_t yval = temp*temp;
    multCount++;
    mimc_constr_wit.add_witness(yval);

    temp = temp*yval;
    multCount++;
    mimc_constr_wit.add_witness(temp);

    xval = temp;


  }

}

#endif
