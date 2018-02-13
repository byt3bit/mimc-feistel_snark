#ifndef VARIABLE_CPP_
#define VARIABLE_CPP_


template<typename field_t>
linear_term< field_t > ::linear_term(const linear_term< field_t >  &other) 
{

  this->coef = other.coef;
  this->index = other.index;
}

template<typename field_t>
linear_term< field_t > & linear_combination< field_t > ::operator[](const size_t i)
{
  if(i > lin_terms.size()) {
    std::cout<<"Error: index out of bound\n";
    return lin_terms[0];
  }
  return lin_terms[i];
}


template<typename field_t>
linear_term< field_t >  linear_term< field_t > ::operator*(const field_t field_coef) 
{
  return linear_term< field_t > (this->index, field_coef * this->coef);
}


template<typename field_t> 
linear_combination< field_t > operator+( const linear_term< field_t > lt1, const linear_term< field_t > lt2 )
{

  linear_combination< field_t > lc;

  if(lt1.index == lt2.index) {

    linear_term< field_t > lt;
    lt.set(lt1.index, lt1.coef + lt2.coef);
    lc.add_term(lt);
    
    
  } else {

    lc.add_term(lt1);
    lc.add_term(lt2);
    
  }

  return lc;
  
}

template< typename field_t > 
linear_combination< field_t > operator+( const linear_combination< field_t > &lc1, const linear_term< field_t > lt2 ) {

  linear_combination< field_t > lc(lc1);

  lc.add_term(lt2);  

  return lc;

}



template<typename field_t>
linear_combination< field_t > ::linear_combination() {}

template< typename field_t >
linear_combination< field_t >::linear_combination(linear_term< field_t > lt )
{
  this->lin_terms.emplace_back(lt);
}


template<typename field_t>
linear_combination< field_t > ::linear_combination(const std::vector<linear_term< field_t >  > &lts)
{
  for(auto t: lts) 
    this->lin_terms.emplace_back(t);
}

template<typename field_t>
linear_combination< field_t > ::linear_combination(const linear_combination< field_t >  &other) {
  
  for(auto t: other.lin_terms) {
    (this->lin_terms).emplace_back(t);
  }
}

template< typename field_t >
linear_combination< field_t > ::linear_combination(std::initializer_list< field_t > lts) {

  for(auto t: lts)
    this->lin_terms.emplace_back(t);
  //this->lin_terms.insert(this->lin_terms.end(), lts.begin(), lts.end());
  
}



template<typename field_t>
linear_combination< field_t >  linear_combination< field_t > ::operator+(const linear_combination< field_t >  &other_lc)
{

  linear_combination< field_t >  result;

  auto *it = this->lin_terms.begin();
  auto *it1 = other_lc.lin_terms.begin();

  while(it != this->lin_terms.end() && it1 != other_lc.lin_terms.end()) {

    if(it->index < it1->index) {
      result.lin_terms.emplace_back(*it);
      ++it;
    }
    else if(it->index > it1->index) {
      result.lin_terms.emplace_back(*it1);
      ++it1;
    }
    else {
      result.lin_terms.emplace_back(linear_term< field_t > (it->index, it->coef+it1->coef));
      ++it;
      ++it1;
    }
  }
  
  while(it != lin_terms.end() ) {
      result.add_term(*it);
      it++;
  }

  while(it1 != lin_terms.end() ) {
      result.add_term(*it1);
      it1++;
  }
  return result;
}


template<typename field_t>
void linear_combination< field_t > ::add_term(const linear_term< field_t >  &lterm) 
{
  this->lin_terms.emplace_back(lterm);
}

template < typename field_t >
void linear_combination < field_t > :: reset(linear_combination < field_t > &other)
{

  this->lin_terms.clear();
  for( auto lt: other)
    lin_terms.emplace_back(lt);
}
  

 

/**newly added for the lowmc function */

// template<typename field_t>
// void linear_combination< field_t > ::append(const svar_array< field_t >  sva)
// {

//   for(auto sv: sva.svars) {
    
//     if(sv.index == 0) { // x_0 --> constant term in the linear combination
//       if(sv.is_assigned) { 
// 	linear_term< field_t >  lt(sv.index, sv.value());
// 	(this->lin_terms).emplace_back(lt);
//       }
//       else {
// 	linear_term< field_t >  lt(sv.index, ONE< field_t > ());
// 	(this->lin_terms).emplace_back(lt);
//       }
//     }
//     else {  
//       linear_term< field_t >  lt(sv.index, ONE< field_t > ());   
//       (this->lin_terms).emplace_back(lt);
//     }
//   }

// }


template<typename field_t>
typename std::vector<linear_term< field_t >  >::const_iterator linear_combination< field_t > ::begin() const
{
  return lin_terms.begin();
}

template<typename field_t>
typename std::vector<linear_term< field_t >  >::const_iterator linear_combination< field_t > ::end() const
{
  return lin_terms.end();
}



template<typename field_t>
using variable_assign = std::vector< field_t > ;

template<typename field_t>
field_t linear_combination< field_t > ::evaluate(std::vector< field_t >  &var_assignment) 
{    
  
  field_t sum;
  for(auto t : lin_terms)     
    sum += (t.index == 0 ? (field_t)1 : t.coef)*var_assignment[t.index-1]; 

  return sum;
}


#endif

