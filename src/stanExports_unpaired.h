// Generated by rstantools.  Do not edit by hand.

/*
    This file is part of rstanarm.
    Copyright (C) 2015, 2016 Trustees of Columbia University
    
    rstanarm is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rstanarm is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with rstanarm.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.26.1-4-gd72b68b7-dirty
#include <stan/model/model_header.hpp>
namespace model_unpaired_namespace {
inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math;
using stan::math::pow; 
stan::math::profile_map profiles__;
static int current_statement__= 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in 'unpaired', line 10, column 2 to column 22)",
                                                      " (in 'unpaired', line 11, column 2 to column 19)",
                                                      " (in 'unpaired', line 12, column 2 to column 30)",
                                                      " (in 'unpaired', line 13, column 2 to column 24)",
                                                      " (in 'unpaired', line 14, column 2 to column 24)",
                                                      " (in 'unpaired', line 17, column 2 to column 19)",
                                                      " (in 'unpaired', line 18, column 2 to column 19)",
                                                      " (in 'unpaired', line 20, column 4 to column 32)",
                                                      " (in 'unpaired', line 19, column 17 to line 21, column 3)",
                                                      " (in 'unpaired', line 19, column 2 to line 21, column 3)",
                                                      " (in 'unpaired', line 23, column 4 to column 39)",
                                                      " (in 'unpaired', line 22, column 17 to line 24, column 3)",
                                                      " (in 'unpaired', line 22, column 2 to line 24, column 3)",
                                                      " (in 'unpaired', line 27, column 2 to column 22)",
                                                      " (in 'unpaired', line 28, column 2 to column 23)",
                                                      " (in 'unpaired', line 29, column 2 to column 20)",
                                                      " (in 'unpaired', line 30, column 2 to column 31)",
                                                      " (in 'unpaired', line 31, column 2 to column 31)",
                                                      " (in 'unpaired', line 32, column 2 to column 31)",
                                                      " (in 'unpaired', line 33, column 2 to column 31)",
                                                      " (in 'unpaired', line 2, column 2 to column 9)",
                                                      " (in 'unpaired', line 3, column 2 to column 9)",
                                                      " (in 'unpaired', line 4, column 16 to column 18)",
                                                      " (in 'unpaired', line 4, column 2 to column 20)",
                                                      " (in 'unpaired', line 5, column 16 to column 18)",
                                                      " (in 'unpaired', line 5, column 2 to column 20)",
                                                      " (in 'unpaired', line 6, column 13 to column 15)",
                                                      " (in 'unpaired', line 6, column 2 to column 17)",
                                                      " (in 'unpaired', line 7, column 12 to column 14)",
                                                      " (in 'unpaired', line 7, column 2 to column 16)",
                                                      " (in 'unpaired', line 13, column 20 to column 22)",
                                                      " (in 'unpaired', line 14, column 20 to column 22)",
                                                      " (in 'unpaired', line 17, column 15 to column 17)",
                                                      " (in 'unpaired', line 18, column 15 to column 17)"};
#include <stan_meta_header.hpp>
class model_unpaired final : public model_base_crtp<model_unpaired> {
private:
  int Ja;
  int Jb;
  std::vector<int> ystararaw;
  std::vector<int> ystarbraw;
  std::vector<double> fpost;
  std::vector<double> fpre;
 
public:
  ~model_unpaired() { }
  
  inline std::string model_name() const final { return "model_unpaired"; }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.26.1-4-gd72b68b7-dirty", "stancflags = "};
  }
  
  
  model_unpaired(stan::io::var_context& context__,
                 unsigned int random_seed__ = 0,
                 std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "model_unpaired_namespace::model_unpaired";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 21;
      context__.validate_dims("data initialization","Ja","int",
          context__.to_vec());
      Ja = std::numeric_limits<int>::min();
      
      current_statement__ = 21;
      Ja = context__.vals_i("Ja")[(1 - 1)];
      current_statement__ = 22;
      context__.validate_dims("data initialization","Jb","int",
          context__.to_vec());
      Jb = std::numeric_limits<int>::min();
      
      current_statement__ = 22;
      Jb = context__.vals_i("Jb")[(1 - 1)];
      current_statement__ = 23;
      validate_non_negative_index("ystararaw", "Ja", Ja);
      current_statement__ = 24;
      context__.validate_dims("data initialization","ystararaw","int",
          context__.to_vec(Ja));
      ystararaw = std::vector<int>(Ja, std::numeric_limits<int>::min());
      
      current_statement__ = 24;
      assign(ystararaw, nil_index_list(), context__.vals_i("ystararaw"),
        "assigning variable ystararaw");
      current_statement__ = 25;
      validate_non_negative_index("ystarbraw", "Jb", Jb);
      current_statement__ = 26;
      context__.validate_dims("data initialization","ystarbraw","int",
          context__.to_vec(Jb));
      ystarbraw = std::vector<int>(Jb, std::numeric_limits<int>::min());
      
      current_statement__ = 26;
      assign(ystarbraw, nil_index_list(), context__.vals_i("ystarbraw"),
        "assigning variable ystarbraw");
      current_statement__ = 27;
      validate_non_negative_index("fpost", "Ja", Ja);
      current_statement__ = 28;
      context__.validate_dims("data initialization","fpost","double",
          context__.to_vec(Ja));
      fpost = std::vector<double>(Ja, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 28;
      assign(fpost, nil_index_list(), context__.vals_r("fpost"),
        "assigning variable fpost");
      current_statement__ = 29;
      validate_non_negative_index("fpre", "Jb", Jb);
      current_statement__ = 30;
      context__.validate_dims("data initialization","fpre","double",
          context__.to_vec(Jb));
      fpre = std::vector<double>(Jb, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 30;
      assign(fpre, nil_index_list(), context__.vals_r("fpre"),
        "assigning variable fpre");
      current_statement__ = 31;
      validate_non_negative_index("mub", "Jb", Jb);
      current_statement__ = 32;
      validate_non_negative_index("mua", "Ja", Ja);
      current_statement__ = 33;
      validate_non_negative_index("lambdaa", "Ja", Ja);
      current_statement__ = 34;
      validate_non_negative_index("lambdab", "Jb", Jb);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      num_params_r__ += 1;
      num_params_r__ += 1;
      num_params_r__ += 1;
      num_params_r__ += Jb;
      num_params_r__ += Ja;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI, stan::require_vector_like_t<VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "model_unpaired_namespace::log_prob";
(void) function__;  // suppress unused var warning
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      local_scalar_t__ kappa;
      kappa = DUMMY_VAR__;
      
      current_statement__ = 1;
      kappa = in__.scalar();
      current_statement__ = 1;
      if (jacobian__) {
        current_statement__ = 1;
        kappa = stan::math::lb_constrain(kappa, 0, lp__);
      } else {
        current_statement__ = 1;
        kappa = stan::math::lb_constrain(kappa, 0);
      }
      local_scalar_t__ mu;
      mu = DUMMY_VAR__;
      
      current_statement__ = 2;
      mu = in__.scalar();
      current_statement__ = 2;
      if (jacobian__) {
        current_statement__ = 2;
        mu = stan::math::lb_constrain(mu, 0, lp__);
      } else {
        current_statement__ = 2;
        mu = stan::math::lb_constrain(mu, 0);
      }
      local_scalar_t__ delta;
      delta = DUMMY_VAR__;
      
      current_statement__ = 3;
      delta = in__.scalar();
      current_statement__ = 3;
      if (jacobian__) {
        current_statement__ = 3;
        delta = stan::math::lub_constrain(delta, 0, 1, lp__);
      } else {
        current_statement__ = 3;
        delta = stan::math::lub_constrain(delta, 0, 1);
      }
      std::vector<local_scalar_t__> mub;
      mub = std::vector<local_scalar_t__>(Jb, DUMMY_VAR__);
      
      current_statement__ = 4;
      for (int sym1__ = 1; sym1__ <= Jb; ++sym1__) {
        current_statement__ = 4;
        assign(mub, cons_list(index_uni(sym1__), nil_index_list()),
          in__.scalar(), "assigning variable mub");}
      current_statement__ = 4;
      for (int sym1__ = 1; sym1__ <= Jb; ++sym1__) {
        current_statement__ = 4;
        if (jacobian__) {
          current_statement__ = 4;
          assign(mub, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(mub[(sym1__ - 1)], 0, lp__),
            "assigning variable mub");
        } else {
          current_statement__ = 4;
          assign(mub, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(mub[(sym1__ - 1)], 0),
            "assigning variable mub");
        }}
      std::vector<local_scalar_t__> mua;
      mua = std::vector<local_scalar_t__>(Ja, DUMMY_VAR__);
      
      current_statement__ = 5;
      for (int sym1__ = 1; sym1__ <= Ja; ++sym1__) {
        current_statement__ = 5;
        assign(mua, cons_list(index_uni(sym1__), nil_index_list()),
          in__.scalar(), "assigning variable mua");}
      current_statement__ = 5;
      for (int sym1__ = 1; sym1__ <= Ja; ++sym1__) {
        current_statement__ = 5;
        if (jacobian__) {
          current_statement__ = 5;
          assign(mua, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(mua[(sym1__ - 1)], 0, lp__),
            "assigning variable mua");
        } else {
          current_statement__ = 5;
          assign(mua, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(mua[(sym1__ - 1)], 0),
            "assigning variable mua");
        }}
      std::vector<local_scalar_t__> lambdaa;
      lambdaa = std::vector<local_scalar_t__>(Ja, DUMMY_VAR__);
      
      std::vector<local_scalar_t__> lambdab;
      lambdab = std::vector<local_scalar_t__>(Jb, DUMMY_VAR__);
      
      current_statement__ = 10;
      for (int i = 1; i <= Jb; ++i) {
        current_statement__ = 8;
        assign(lambdab, cons_list(index_uni(i), nil_index_list()),
          (mub[(i - 1)] / fpre[(i - 1)]), "assigning variable lambdab");}
      current_statement__ = 13;
      for (int i = 1; i <= Ja; ++i) {
        current_statement__ = 11;
        assign(lambdaa, cons_list(index_uni(i), nil_index_list()),
          ((delta * mua[(i - 1)]) / fpost[(i - 1)]),
          "assigning variable lambdaa");}
      {
        current_statement__ = 14;
        lp_accum__.add(gamma_lpdf<propto__>(mu, 1, 0.001));
        current_statement__ = 15;
        lp_accum__.add(gamma_lpdf<propto__>(kappa, 1, 0.7));
        current_statement__ = 16;
        lp_accum__.add(beta_lpdf<propto__>(delta, 1, 1));
        current_statement__ = 17;
        lp_accum__.add(gamma_lpdf<propto__>(mub, kappa, (kappa / mu)));
        current_statement__ = 18;
        lp_accum__.add(gamma_lpdf<propto__>(mua, kappa, (kappa / mu)));
        current_statement__ = 19;
        lp_accum__.add(poisson_lpmf<propto__>(ystarbraw, lambdab));
        current_statement__ = 20;
        lp_accum__.add(poisson_lpmf<propto__>(ystararaw, lambdaa));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr>
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "model_unpaired_namespace::write_array";
(void) function__;  // suppress unused var warning
    (void) function__;  // suppress unused var warning
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      double kappa;
      kappa = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      kappa = in__.scalar();
      current_statement__ = 1;
      kappa = stan::math::lb_constrain(kappa, 0);
      double mu;
      mu = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      mu = in__.scalar();
      current_statement__ = 2;
      mu = stan::math::lb_constrain(mu, 0);
      double delta;
      delta = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      delta = in__.scalar();
      current_statement__ = 3;
      delta = stan::math::lub_constrain(delta, 0, 1);
      std::vector<double> mub;
      mub = std::vector<double>(Jb, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 4;
      for (int sym1__ = 1; sym1__ <= Jb; ++sym1__) {
        current_statement__ = 4;
        assign(mub, cons_list(index_uni(sym1__), nil_index_list()),
          in__.scalar(), "assigning variable mub");}
      current_statement__ = 4;
      for (int sym1__ = 1; sym1__ <= Jb; ++sym1__) {
        current_statement__ = 4;
        assign(mub, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_constrain(mub[(sym1__ - 1)], 0),
          "assigning variable mub");}
      std::vector<double> mua;
      mua = std::vector<double>(Ja, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 5;
      for (int sym1__ = 1; sym1__ <= Ja; ++sym1__) {
        current_statement__ = 5;
        assign(mua, cons_list(index_uni(sym1__), nil_index_list()),
          in__.scalar(), "assigning variable mua");}
      current_statement__ = 5;
      for (int sym1__ = 1; sym1__ <= Ja; ++sym1__) {
        current_statement__ = 5;
        assign(mua, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_constrain(mua[(sym1__ - 1)], 0),
          "assigning variable mua");}
      std::vector<double> lambdaa;
      lambdaa = std::vector<double>(Ja, std::numeric_limits<double>::quiet_NaN());
      
      std::vector<double> lambdab;
      lambdab = std::vector<double>(Jb, std::numeric_limits<double>::quiet_NaN());
      
      vars__.emplace_back(kappa);
      vars__.emplace_back(mu);
      vars__.emplace_back(delta);
      for (int sym1__ = 1; sym1__ <= Jb; ++sym1__) {
        vars__.emplace_back(mub[(sym1__ - 1)]);}
      for (int sym1__ = 1; sym1__ <= Ja; ++sym1__) {
        vars__.emplace_back(mua[(sym1__ - 1)]);}
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      current_statement__ = 10;
      for (int i = 1; i <= Jb; ++i) {
        current_statement__ = 8;
        assign(lambdab, cons_list(index_uni(i), nil_index_list()),
          (mub[(i - 1)] / fpre[(i - 1)]), "assigning variable lambdab");}
      current_statement__ = 13;
      for (int i = 1; i <= Ja; ++i) {
        current_statement__ = 11;
        assign(lambdaa, cons_list(index_uni(i), nil_index_list()),
          ((delta * mua[(i - 1)]) / fpost[(i - 1)]),
          "assigning variable lambdaa");}
      if (emit_transformed_parameters__) {
        for (int sym1__ = 1; sym1__ <= Ja; ++sym1__) {
          vars__.emplace_back(lambdaa[(sym1__ - 1)]);}
        for (int sym1__ = 1; sym1__ <= Jb; ++sym1__) {
          vars__.emplace_back(lambdab[(sym1__ - 1)]);}
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, stan::require_std_vector_t<VecVar>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void transform_inits_impl(const stan::io::var_context& context__,
                                   VecI& params_i__, VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.clear();
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      double kappa;
      kappa = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      kappa = context__.vals_r("kappa")[(1 - 1)];
      double kappa_free__;
      kappa_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      kappa_free__ = stan::math::lb_free(kappa, 0);
      double mu;
      mu = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      mu = context__.vals_r("mu")[(1 - 1)];
      double mu_free__;
      mu_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      mu_free__ = stan::math::lb_free(mu, 0);
      double delta;
      delta = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      delta = context__.vals_r("delta")[(1 - 1)];
      double delta_free__;
      delta_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      delta_free__ = stan::math::lub_free(delta, 0, 1);
      std::vector<double> mub;
      mub = std::vector<double>(Jb, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 4;
      assign(mub, nil_index_list(), context__.vals_r("mub"),
        "assigning variable mub");
      std::vector<double> mub_free__;
      mub_free__ = std::vector<double>(Jb, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 4;
      for (int sym1__ = 1; sym1__ <= Jb; ++sym1__) {
        current_statement__ = 4;
        assign(mub_free__, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_free(mub[(sym1__ - 1)], 0),
          "assigning variable mub_free__");}
      std::vector<double> mua;
      mua = std::vector<double>(Ja, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 5;
      assign(mua, nil_index_list(), context__.vals_r("mua"),
        "assigning variable mua");
      std::vector<double> mua_free__;
      mua_free__ = std::vector<double>(Ja, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 5;
      for (int sym1__ = 1; sym1__ <= Ja; ++sym1__) {
        current_statement__ = 5;
        assign(mua_free__, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_free(mua[(sym1__ - 1)], 0),
          "assigning variable mua_free__");}
      vars__.emplace_back(kappa_free__);
      vars__.emplace_back(mu_free__);
      vars__.emplace_back(delta_free__);
      for (int sym1__ = 1; sym1__ <= Jb; ++sym1__) {
        vars__.emplace_back(mub_free__[(sym1__ - 1)]);}
      for (int sym1__ = 1; sym1__ <= Ja; ++sym1__) {
        vars__.emplace_back(mua_free__[(sym1__ - 1)]);}
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__.clear();
    names__.emplace_back("kappa");
    names__.emplace_back("mu");
    names__.emplace_back("delta");
    names__.emplace_back("mub");
    names__.emplace_back("mua");
    names__.emplace_back("lambdaa");
    names__.emplace_back("lambdab");
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.clear();
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(Jb)});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(Ja)});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(Ja)});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(Jb)});
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "kappa");
    param_names__.emplace_back(std::string() + "mu");
    param_names__.emplace_back(std::string() + "delta");
    for (int sym1__ = 1; sym1__ <= Jb; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "mub" + '.' + std::to_string(sym1__));
      }}
    for (int sym1__ = 1; sym1__ <= Ja; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "mua" + '.' + std::to_string(sym1__));
      }}
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= Ja; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "lambdaa" + '.' + std::to_string(sym1__));
        }}
      for (int sym1__ = 1; sym1__ <= Jb; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "lambdab" + '.' + std::to_string(sym1__));
        }}
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "kappa");
    param_names__.emplace_back(std::string() + "mu");
    param_names__.emplace_back(std::string() + "delta");
    for (int sym1__ = 1; sym1__ <= Jb; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "mub" + '.' + std::to_string(sym1__));
      }}
    for (int sym1__ = 1; sym1__ <= Ja; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "mua" + '.' + std::to_string(sym1__));
      }}
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= Ja; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "lambdaa" + '.' + std::to_string(sym1__));
        }}
      for (int sym1__ = 1; sym1__ <= Jb; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "lambdab" + '.' + std::to_string(sym1__));
        }}
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"kappa\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"delta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"mub\",\"type\":{\"name\":\"array\",\"length\":" << Jb << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"mua\",\"type\":{\"name\":\"array\",\"length\":" << Ja << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"lambdaa\",\"type\":{\"name\":\"array\",\"length\":" << Ja << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"lambdab\",\"type\":{\"name\":\"array\",\"length\":" << Jb << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"kappa\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"delta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"mub\",\"type\":{\"name\":\"array\",\"length\":" << Jb << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"mua\",\"type\":{\"name\":\"array\",\"length\":" << Ja << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"lambdaa\",\"type\":{\"name\":\"array\",\"length\":" << Ja << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"lambdab\",\"type\":{\"name\":\"array\",\"length\":" << Jb << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      std::vector<double> vars_vec(vars.size());
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i) {
        vars.coeffRef(i) = vars_vec[i];
      }
    }
    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }
    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
  
    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits_impl(context, params_i, params_r_vec, pstream);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i) {
        params_r.coeffRef(i) = params_r_vec[i];
      }
    }
    inline void transform_inits(const stan::io::var_context& context,
                                std::vector<int>& params_i,
                                std::vector<double>& vars,
                                std::ostream* pstream = nullptr) const final {
      transform_inits_impl(context, params_i, vars, pstream);
    }        
};
}
using stan_model = model_unpaired_namespace::model_unpaired;
#ifndef USING_R
// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_unpaired_namespace::profiles__;
}
#endif
#endif
