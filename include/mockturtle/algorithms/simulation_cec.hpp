/* mockturtle: C++ logic network library

 * Copyright (C) 2018-2021  EPFL

 *

 * Permission is hereby granted, free of charge, to any person

 * obtaining a copy of this software and associated documentation

 * files (the "Software"), to deal in the Software without

 * restriction, including without limitation the rights to use,

 * copy, modify, merge, publish, distribute, sublicense, and/or sell

 * copies of the Software, and to permit persons to whom the

 * Software is furnished to do so, subject to the following

 * conditions:

 *

 * The above copyright notice and this permission notice shall be

 * included in all copies or substantial portions of the Software.

 *

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,

 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES

 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND

 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT

 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,

 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING

 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR

 * OTHER DEALINGS IN THE SOFTWARE.

 */



/*!

  \file simulation_cec.hpp

  \brief Simulation-based CEC



  EPFL CS-472 2021 Final Project Option 2

*/



#pragma once



#include <kitty/constructors.hpp>

#include <kitty/dynamic_truth_table.hpp>

#include <kitty/operations.hpp>



#include "../utils/node_map.hpp"

#include "miter.hpp"

#include"math.h" 

#include "simulation.hpp"

#include <iostream>
#include <cmath>

namespace mockturtle


{



/* Statistics to be reported */

struct simulation_cec_stats

{

  /*! \brief Split variable (simulation size). */

  uint32_t split_var{ 0 };



  /*! \brief Number of simulation rounds. */

  uint32_t rounds{ 0 };

};



namespace detail

{



template<class Ntk>

class simulation_cec_impl

{

public:

  using pattern_t = unordered_node_map<kitty::dynamic_truth_table, Ntk>;

  using node = typename Ntk::node;

  using signal = typename Ntk::signal;



public:

  explicit simulation_cec_impl( Ntk& ntk, simulation_cec_stats& st )

      : _ntk(ntk),

        _st(st)

  {

  }



  bool run()

  {

  _st.split_var = computing_splitting_var(_ntk.num_pis());
  _st.rounds = compute_round( _ntk.num_pis(), _st.split_var);


  
  pattern_t pattern(_ntk);
  default_simulator <kitty::dynamic_truth_table> sim(_st.split_var);
   init_patterns (pattern);
   simulate_nodes(_ntk, pattern ,sim);
   if (check_patterns(pattern)==false){
     return false;
   }
  for (uint32_t i=1; i<=_st.rounds-1 ; i++){
    remove_pattern( pattern );
    update_patterns(  pattern,  i);
    simulate_nodes(_ntk, pattern ,sim);
    if (check_patterns(pattern)==false){
      return false;
    }
  }
  return true;

  }



private:

   uint32_t computing_splitting_var(uint32_t n){
        float k;
         k=(log((pow(2,29)/_ntk._storage->nodes.size()) - 32)/log(2))+3;
          if (n<7){
           _st.split_var = n;
         }
          else {
            if (n<=k){_st.split_var=n;}
            else {_st.split_var=k;}
         }
         return _st.split_var;
        
         }
 uint32_t compute_round(uint32_t n, uint32_t split_var){

        uint32_t rounds;
        rounds = pow(2,n-split_var);
        return rounds;
      }

   void init_patterns( pattern_t& pattern){
    _ntk.foreach_pi( [&]( auto const& n, auto i ) { 
          kitty::dynamic_truth_table tt(_st.split_var);
          if (i<_st.split_var){
            kitty::create_nth_var(tt,i);}

            pattern[n]= tt;
          
  });
  }
    bool check_patterns( pattern_t& pattern){
      bool equiv;
      equiv = true;
      _ntk.foreach_po([&]( auto const& j ){
        if (_ntk.is_complemented(j)){
          if (!is_const0(~pattern[j])){ 
            equiv=false;
          }
        }
          else {
            if (!is_const0(pattern[j])) {
              equiv=false;
            }
          }  
        });
      return equiv;
    }
    
  

    void update_patterns (  pattern_t& pattern,  uint32_t i){
      uint32_t p = i;
      _ntk.foreach_pi( [&]( auto const& n, auto j ) {
          if (j>= _st.split_var){
            if (p%2 ==1 ){
              if (is_const0(pattern[n])){
                pattern[n]=~pattern[n];
              }
            }
            else{
              if (!is_const0(pattern[n])){
                pattern[n]=~pattern[n];
              }
            }
            p=p/2;
          }
        });
    }
    
    void remove_pattern(pattern_t& pattern){
      _ntk.foreach_gate([&](auto const& i){
        pattern.erase(i);
      });
    }


private:

  Ntk& _ntk;

  simulation_cec_stats& _st;

  /* you can add other attributes here */

};



} // namespace detail



/* Entry point for users to call */



/*! \brief Simulation-based CEC.

 *

 * This function implements a simulation-based combinational equivalence checker.

 * The implementation creates a miter network and run several rounds of simulation

 * to verify the functional equivalence. For memory and speed reasons this approach

 * is limited up to 40 input networks. It returns an optional which is `nullopt`,

 * if the network has more than 40 inputs.

 */

template<class Ntk>

std::optional<bool> simulation_cec( Ntk const& ntk1, Ntk const& ntk2, simulation_cec_stats* pst = nullptr )

{

  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );

  static_assert( has_num_pis_v<Ntk>, "Ntk does not implement the num_pis method" );

  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );

  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );

  static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );

  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );

  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );

  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );



  simulation_cec_stats st;



  bool result = false;



  if ( ntk1.num_pis() > 40 )

    return std::nullopt;



  auto ntk_miter = miter<Ntk>( ntk1, ntk2 );



  if ( ntk_miter.has_value() )

  {

    detail::simulation_cec_impl p( *ntk_miter, st );

    result = p.run();

  }



  if ( pst )

    *pst = st;



  return result;

}



} // namespace mockturtle