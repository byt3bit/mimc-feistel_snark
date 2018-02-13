#ifndef MIMC_FSTL_2BR_SNARK_CPP
#define MIMC_FSTL_2BR_SNARK_CPP

/*

This function computes all the rank-1 constraints for the MiMC 
2 branch FN                           
+--------------------------------------------------------------------------+
*/

template<typename field_t> 
void mimcfstl_2br_snark<field_t>::generate_r1_constraint()
{


  index_t lindex = 1;
  index_t rindex = lindex + 1;
  index_t var_index = rindex;

  var_index++;
  for(int i = 0;i < numofRound;i++) {

    linear_term< field_t > xc(0, roundConst[i]);
    linear_term< field_t > x1(lindex, (field_t)(1));
    linear_term< field_t > x2(rindex, (field_t)(1));
    
    linear_term< field_t > x3(var_index, (field_t)(1));
    var_index++;
    linear_combination< field_t > A(xc + x1);
    linear_combination < field_t > B(A);
    linear_combination < field_t > C(x3);

    constraint < field_t > constr(A, B, C);
    mimcfstl2br_constr_wit.add_constraint(constr);

    A.reset(C);
    x3.set(var_index, (field_t) 1);
    x2.set(rindex, (field_t) (-1) );
    rindex = lindex;
    lindex = var_index;
   
    var_index++;
    
    C.clear();
    C.add_term(x3);
    C.add_term(x2);
    constr.reset_constraint(A, B, C);
    
    mimcfstl2br_constr_wit.add_constraint(constr);

  }

  //mimcfstl2br_constr_wit.print_constraints();

}


template<typename field_t> 
void mimcfstl_2br_snark<field_t>::generate_witness(field_t leftinput,
						   field_t rightinput)
{
  
  index_t var_index = 0;
  mimcfstl2br_constr_wit.add_witness({leftinput, rightinput});
  var_index++;
  field_t lval = leftinput;
  field_t rval = rightinput;

  for(int i = 0;i < numofRound;i++) {

    field_t temp = (leftinput + roundConst[i]);
    field_t temp1 = temp*temp;
    addCount += 1;
    multCount++;
    
    //field_t temp1 = (leftinput + roundConst[i])*temp + rightinput;
    temp = temp1*temp + rightinput;
    //temp = temp + rightinput; 
    addCount += 1;
    multCount++;

    mimcfstl2br_constr_wit.add_witness({temp, temp1});

    rightinput = leftinput;
    leftinput = temp1;

  }

}


#endif
