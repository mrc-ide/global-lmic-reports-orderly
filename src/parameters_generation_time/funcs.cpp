// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/gamma.hpp>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  // per block setup
  double prop_Mild;
  double gamma_E;
  double gamma_ICase;
  double gamma_IMild;

  Rcpp::NumericVector gens;
  List model_list = misc["model_instance"];
  Function set_user = Rcpp::as<Function>(model_list["set_user"]);
  Function run_model = Rcpp::as<Function>(model_list["run"]);
  int block_num = misc["block"];
  Rcpp::NumericVector ts = misc["denom_values"];
  int n_ts = ts.length();
  double trapezoid_multipler = misc["trapezoid_multipler"];

  if(block_num == 1){
    prop_Mild = 0.26;
    gamma_E = 2/params["dur_E"];
    gamma_ICase = 2/params["dur_ICase"];
    gamma_IMild = 1/params["dur_IMild"];
    gens = data["wild"];
  } else if(block_num == 2){
    prop_Mild = 0.9426399;
    gamma_E = 2/params["dur_E_delta"];
    gamma_ICase = 2/params["dur_ICase_delta"];
    gamma_IMild = 1/params["dur_IMild_delta"];
    gens = data["delta"];
  } else {
    prop_Mild = 0.9426399;
    gamma_E = 2/params["dur_E_omicron"];
    gamma_ICase = 2/params["dur_ICase_omicron"];
    gamma_IMild = 1/params["dur_IMild_omicron"];
    gens = data["omicron"];
  }
  //apply parameters to odin model
  set_user(Named("prop_Mild", prop_Mild),
           Named("gamma_E", gamma_E),
           Named("gamma_ICase", gamma_ICase),
           Named("gamma_IMild", gamma_IMild));
  //run the model only twice to improve performance
  //must be some way to combine these easily
  Rcpp::NumericMatrix model_run = run_model(Named("t", ts));
  Rcpp::NumericVector denom_infectiousness = model_run(_, 6);
  //get denominator via trapeziod
  double denom = 0;
  for(R_xlen_t i = 1; i < n_ts - 1; i++) {
    double new_value = denom_infectiousness[i];
    if(new_value < 0.0001){
      new_value = 0.0001;
    }
    denom += new_value;
  }
  denom += denom + denom_infectiousness[0] + denom_infectiousness[n_ts - 1];
  double log_denom = log(denom) + trapezoid_multipler;
  //get infectiousness
  gens.push_front(0);
  Rcpp::NumericMatrix model_run_2 = run_model(Named("t", gens));
  Rcpp::NumericVector gen_infectiousness = model_run_2(_, 6);
  gen_infectiousness.erase(0);
  //calculate likelihood
  double loglikelihood = 0;
  for(int i = 0; i < gen_infectiousness.length(); i++){
    //set a minimum so we don't get NAN values from negatives
    double log_value = gen_infectiousness[i];
    if(log_value < 0.0001){
      log_value = 0.0001;
    }
    loglikelihood += log(log_value);
  }
  return Rcpp::wrap(loglikelihood - log_denom);
}

// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  //parameters
  float dur_E = params["dur_E"];
  float dur_ICase = params["dur_ICase"];
  float dur_IMild = params["dur_IMild"];
  float dur_E_delta = params["dur_E_delta"];
  float dur_ICase_delta = params["dur_ICase_delta"];
  float dur_IMild_delta = params["dur_IMild_delta"];
  float dur_E_omicron = params["dur_E_omicron"];
  float dur_ICase_omicron = params["dur_ICase_omicron"];
  float dur_IMild_omicron = params["dur_IMild_omicron"];
  //distributions
  boost::math::gamma_distribution<> dur_E_dist(2, 2.3);
  boost::math::gamma_distribution<> dur_E_ICase_dist(2, 2.25);
  boost::math::gamma_distribution<> dur_E_IMild_dist(1, 2.1);
  boost::math::normal_distribution<> changes_dist(0, 1);
  //calculate prior, logpdf not yet supported?
  double logprior = log(pdf(dur_E_dist, dur_E)) +
    log(pdf(dur_E_ICase_dist, dur_ICase)) +
    log(pdf(dur_E_IMild_dist, dur_IMild)) +
    log(pdf(changes_dist, dur_E_delta - dur_E)) +
    log(pdf(changes_dist, dur_ICase_delta - dur_ICase)) +
    log(pdf(changes_dist, dur_IMild_delta - dur_IMild)) +
    log(pdf(changes_dist, dur_E_omicron - dur_E)) +
    log(pdf(changes_dist, dur_ICase_omicron - dur_ICase)) +
    log(pdf(changes_dist, dur_IMild_omicron - dur_IMild));
  return Rcpp::wrap(logprior);
}


// NOTE: Do not edit this function name
// [[Rcpp::export]]
SEXP create_xptr(std::string function_name) {
  typedef SEXP (*funcPtr_likelihood)(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc);
  typedef SEXP (*funcPtr_prior)(Rcpp::NumericVector params, Rcpp::List misc);

  // NOTE: If your loglikelihood function is not called "loglike" please edit:
  if (function_name == "loglike"){
    return(Rcpp::XPtr<funcPtr_likelihood>(new funcPtr_likelihood(&loglike)));
  }
  // NOTE: If your logprior function is not called "logprior" please edit:
  if (function_name == "logprior"){
    return(Rcpp::XPtr<funcPtr_prior>(new funcPtr_prior(&logprior)));
  }

  stop("cpp function %i not found", function_name);
}
