//
// vegetation predictive process
//

#define EIGEN_DONT_PARALLELIZE

#include <stan/version.hpp>
#include <stan/io/cmd_line.hpp>
#include <stan/io/dump.hpp>
#include <stan/io/mcmc_writer.hpp>
#include <stan/io/writer.hpp>
#include <stan/io/reader.hpp>

#include <stan/model/prob_grad.hpp>
#include <stan/math/prim.hpp>

#include <stan/services/io.hpp>
#include <stan/services/mcmc.hpp>

#include <string>
#include <vector>

#include <unsupported/Eigen/KroneckerProduct>

#include "timer.hpp"

namespace veg_pp {

  using namespace std;
  using namespace stan::math;
  using namespace stan::prob;
  using namespace stan::io;

  typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
  typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
  typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

  class veg_pp_model : public stan::model::prob_grad {
  private:
    int K;
    int N;
    int N_knots;
    vector< vector<int> > y;
    matrix_d d_knots;
    matrix_d d_inter;
    int N_p;
    double P;
    vector_d ones;
    vector_d zeros;
    int W;
    matrix_d eye;
    //matrix_d M;
    Eigen::HouseholderQR<vector_d> qr;
    matrix_d fatQ;
    vector_d Q;
    row_vector_d QT;

  public:
    veg_pp_model(stan::io::var_context& context__, std::ostream* pstream__ = 0) : prob_grad::prob_grad(0) {

      static const char* function__ = "veg_pp::veg_pp_model(%1%)";
      size_t pos__;
      std::vector<int> vals_i__;
      std::vector<double> vals_r__;

      context__.validate_dims("data initialization", "K", "int", context__.to_vec());
      K = int(0);
      vals_i__ = context__.vals_i("K");
      pos__ = 0;
      K = vals_i__[pos__++];

      context__.validate_dims("data initialization", "N", "int", context__.to_vec());
      N = int(0);
      vals_i__ = context__.vals_i("N");
      pos__ = 0;
      N = vals_i__[pos__++];

      context__.validate_dims("data initialization", "N_knots", "int", context__.to_vec());
      N_knots = int(0);
      vals_i__ = context__.vals_i("N_knots");
      pos__ = 0;
      N_knots = vals_i__[pos__++];

      context__.validate_dims("data initialization", "y", "int", context__.to_vec(N,K));
      stan::math::validate_non_negative_index("y", "N", N);
      stan::math::validate_non_negative_index("y", "K", K);
      y = std::vector<std::vector<int> >(N,std::vector<int>(K,int(0)));
      vals_i__ = context__.vals_i("y");
      pos__ = 0;
      size_t y_limit_1__ = K;
      for (size_t i_1__ = 0; i_1__ < y_limit_1__; ++i_1__) {
	size_t y_limit_0__ = N;
	for (size_t i_0__ = 0; i_0__ < y_limit_0__; ++i_0__) {
	  y[i_0__][i_1__] = vals_i__[pos__++];
	}
      }

      context__.validate_dims("data initialization", "d_knots", "matrix_d", context__.to_vec(N_knots,N_knots));
      stan::math::validate_non_negative_index("d_knots", "N_knots", N_knots);
      stan::math::validate_non_negative_index("d_knots", "N_knots", N_knots);
      d_knots = matrix_d(N_knots,N_knots);
      vals_r__ = context__.vals_r("d_knots");
      pos__ = 0;
      size_t d_knots_m_mat_lim__ = N_knots;
      size_t d_knots_n_mat_lim__ = N_knots;
      for (size_t n_mat__ = 0; n_mat__ < d_knots_n_mat_lim__; ++n_mat__) {
	for (size_t m_mat__ = 0; m_mat__ < d_knots_m_mat_lim__; ++m_mat__) {
	  d_knots(m_mat__,n_mat__) = vals_r__[pos__++];
	}
      }

      context__.validate_dims("data initialization", "d_inter", "matrix_d", context__.to_vec(N,N_knots));
      stan::math::validate_non_negative_index("d_inter", "N", N);
      stan::math::validate_non_negative_index("d_inter", "N_knots", N_knots);
      d_inter = matrix_d(N,N_knots);
      vals_r__ = context__.vals_r("d_inter");
      pos__ = 0;
      size_t d_inter_m_mat_lim__ = N;
      size_t d_inter_n_mat_lim__ = N_knots;
      for (size_t n_mat__ = 0; n_mat__ < d_inter_n_mat_lim__; ++n_mat__) {
	for (size_t m_mat__ = 0; m_mat__ < d_inter_m_mat_lim__; ++m_mat__) {
	  d_inter(m_mat__,n_mat__) = vals_r__[pos__++];
	}
      }

      context__.validate_dims("data initialization", "N_p", "int", context__.to_vec());
      N_p = int(0);
      vals_i__ = context__.vals_i("N_p");
      pos__ = 0;
      N_p = vals_i__[pos__++];

      context__.validate_dims("data initialization", "P", "double", context__.to_vec());
      P = double(0);
      vals_r__ = context__.vals_r("P");
      pos__ = 0;
      P = vals_r__[pos__++];

      try {
	check_greater_or_equal(function__,"K",K,0);
      } catch (std::domain_error& e) {
	throw std::domain_error(std::string("Invalid value of K: ") + std::string(e.what()));
      };
      try {
	check_greater_or_equal(function__,"N",N,0);
      } catch (std::domain_error& e) {
	throw std::domain_error(std::string("Invalid value of N: ") + std::string(e.what()));
      };
      try {
	check_greater_or_equal(function__,"N_knots",N_knots,0);
      } catch (std::domain_error& e) {
	throw std::domain_error(std::string("Invalid value of N_knots: ") + std::string(e.what()));
      };
      try {
	check_greater_or_equal(function__,"N_p",N_p,0);
      } catch (std::domain_error& e) {
	throw std::domain_error(std::string("Invalid value of N_p: ") + std::string(e.what()));
      };

      // transformed data
      stan::math::validate_non_negative_index("ones", "N", N);
      ones = vector_d(N);
      stan::math::validate_non_negative_index("zeros", "N_knots", N_knots);
      zeros = vector_d(N_knots);
      W = int(0);

      // eye = matrix_d(N_knots,N_knots).setIdentity();

      // stan::math::validate_non_negative_index("M", "N", N);
      // M = matrix_d(N,N);

      // for (int j = 1; j <= N; ++j) {
      // 	for (int i = 1; i <= N; ++i) {
      // 	  stan::math::assign(get_base1_lhs(M,i,j,"M",1), -(P));
      // 	}
      // 	stan::math::assign(get_base1_lhs(M,j,j,"M",1), (1.0 - P));
      // }
      for (int i = 1; i <= N; ++i) {
	stan::math::assign(get_base1_lhs(ones,i,"ones",1), 1.0);
      }
      for (int i = 1; i <= N_knots; ++i) {
	stan::math::assign(get_base1_lhs(zeros,i,"zeros",1), 0.0);
      }
      //stan::math::assign(K, (K - 1));

      qr = ones.householderQr();
      fatQ = qr.householderQ();
      Q = fatQ.col(0); // thin Q, as returned in R
      QT = Q.transpose();

      // set parameter ranges
      num_params_r__ = 0U;
      param_ranges_i__.clear();
      num_params_r__ += K;
      num_params_r__ += K;
      num_params_r__ += K;
      num_params_r__ += N_knots * K;
      num_params_r__ += N * K;
    }

    ~veg_pp_model() { }

    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__) const {
      stan::io::writer<double> writer__(params_r__,params_i__);
      size_t pos__;
      (void) pos__; // dummy call to supress warning
      std::vector<double> vals_r__;
      std::vector<int> vals_i__;


      if (!(context__.contains_r("eta")))
	throw std::runtime_error("variable eta missing");
      vals_r__ = context__.vals_r("eta");
      pos__ = 0U;
      context__.validate_dims("initialization", "eta", "vector_d", context__.to_vec(K));
      vector_d eta(K);
      for (int j1__ = 0U; j1__ < K; ++j1__)
	eta(j1__) = vals_r__[pos__++];
      try {
	writer__.vector_lub_unconstrain(0.10000000000000001,sqrt(100),eta);
      } catch (std::exception& e) {
	throw std::runtime_error(std::string("Error transforming variable eta: ") + e.what());
      }

      if (!(context__.contains_r("rho")))
	throw std::runtime_error("variable rho missing");
      vals_r__ = context__.vals_r("rho");
      pos__ = 0U;
      context__.validate_dims("initialization", "rho", "vector_d", context__.to_vec(K));
      vector_d rho(K);
      for (int j1__ = 0U; j1__ < K; ++j1__)
	rho(j1__) = vals_r__[pos__++];
      try {
	writer__.vector_lub_unconstrain(1.0e-6,1.0,rho);
      } catch (std::exception& e) {
	throw std::runtime_error(std::string("Error transforming variable rho: ") + e.what());
      }

      if (!(context__.contains_r("mu")))
	throw std::runtime_error("variable mu missing");
      vals_r__ = context__.vals_r("mu");
      pos__ = 0U;
      context__.validate_dims("initialization", "mu", "vector_d", context__.to_vec(K));
      vector_d mu(K);
      for (int j1__ = 0U; j1__ < K; ++j1__)
	mu(j1__) = vals_r__[pos__++];
      try {
	writer__.vector_unconstrain(mu);
      } catch (std::exception& e) {
	throw std::runtime_error(std::string("Error transforming variable mu: ") + e.what());
      }

        if (!(context__.contains_r("alpha")))
            throw std::runtime_error("variable alpha missing");
        vals_r__ = context__.vals_r("alpha");
        pos__ = 0U;
        context__.validate_dims("initialization", "alpha", "vector_d", context__.to_vec(K,N_knots));
        std::vector<vector_d> alpha(K,vector_d(N_knots));
        for (int j1__ = 0U; j1__ < N_knots; ++j1__)
            for (int i0__ = 0U; i0__ < K; ++i0__)
                alpha[i0__](j1__) = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < K; ++i0__)
            try { writer__.vector_unconstrain(alpha[i0__]); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable alpha: ") + e.what()); }

        if (!(context__.contains_r("g")))
            throw std::runtime_error("variable g missing");
        vals_r__ = context__.vals_r("g");
        pos__ = 0U;
        context__.validate_dims("initialization", "g", "vector_d", context__.to_vec(K,N));
        std::vector<vector_d> g(K,vector_d(N));
        for (int j1__ = 0U; j1__ < N; ++j1__)
            for (int i0__ = 0U; i0__ < K; ++i0__)
                g[i0__](j1__) = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < K; ++i0__)
            try { writer__.vector_unconstrain(g[i0__]); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable g: ") + e.what()); }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

   void transform_inits(const stan::io::var_context& context,
                         Eigen::VectorXd& params_r,
                         ostream* output) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    inline double lub_transform(const double x, const double lb, const double ub,
			 double &lja, double &ja, double &dj) const
    {
      double inv_logit_x;
      if (x > 0) {
        double exp_minus_x = exp(-x);
	double exp_minus_x_p1 = exp_minus_x + 1.0;
        inv_logit_x = 1.0 / (1.0 + exp_minus_x);
        lja = log(ub - lb) - x - 2 * log1p(exp_minus_x);
	ja  = (ub - lb) * exp_minus_x / (exp_minus_x_p1 * exp_minus_x_p1);
	dj  = -1.0 + 2.0 * exp_minus_x / (1 + exp_minus_x);
        if ((x < std::numeric_limits<double>::infinity())
            && (inv_logit_x == 1))
            inv_logit_x = 1 - 1e-15;
      } else {
        double exp_x = exp(x);
	double exp_x_p1 = exp_x + 1.0;
        inv_logit_x = 1.0 - 1.0 / (1.0 + exp_x);
        lja = log(ub - lb) + x - 2 * log1p(exp_x);
	ja  = (ub - lb) * exp_x / (exp_x_p1 * exp_x_p1);
	dj  = -1.0 + 2.0 / (exp_x + 1.0);
        if ((x > -std::numeric_limits<double>::infinity())
            && (inv_logit_x== 0))
            inv_logit_x = 1e-15;
      }
      return lb + (ub - lb) * inv_logit_x;
    }

    double normal_log_double(const vector_d& y, const double mu, const double sigma) const
    {
      double lp = 0.0;
      double inv_sigma = 1.0/sigma;
      double log_sigma = log(sigma);

      // #ifndef NDEBUG
      // if (!check_finite("normal_log_double", mu, "Location parameter", &lp))
      //   return lp;
      // if (!check_not_nan("normal_log_double", y, "Random variable", &lp))
      //   return lp;
      // #endif

      for (int n = 0; n < y.size(); n++) {
	const double y_minus_mu_over_sigma = (y[n] - mu) * inv_sigma;
	const double y_minus_mu_over_sigma_squared = y_minus_mu_over_sigma * y_minus_mu_over_sigma;
	lp -= 0.5 * y_minus_mu_over_sigma_squared;
      }
      return lp;
    }

    double normal_log_double(const double y, const double mu, const double sigma) const
    {
      double lp = 0.0;
      double inv_sigma = 1.0/sigma;
      double log_sigma = log(sigma);

      const double y_minus_mu_over_sigma = (y - mu) * inv_sigma;
      const double y_minus_mu_over_sigma_squared = y_minus_mu_over_sigma * y_minus_mu_over_sigma;
      lp -= 0.5 * y_minus_mu_over_sigma_squared;
      lp -= log_sigma;

      return lp;
    }

    double multi_normal_cholesky_log_double(const vector_d& y, const vector_d& mu, const matrix_d& L) const {

      double lp = 0.0;

      // #ifndef NDEBUG
      // if (!check_finite("multi_normal_cholesky_log_double", mu, "Location parameter", &lp))
      //   return lp;
      // if (!check_not_nan("multi_normal_cholesky_log_double(%1%)", y, "Random variable", &lp))
      //   return lp;
      // #endif

      vector_d L_log_diag = L.diagonal().array().log().matrix();
      lp -= sum(L_log_diag);

      vector_d y_minus_mu = y - mu;
      vector_d half = mdivide_left_tri_low(L, y_minus_mu);
      lp -= 0.5 * dot_self(half);

      return lp;
    }

    double multinomial_log_double(const std::vector<int>& ns,
				  const vector_d& theta) const {

      double lp = 0.0;
      // #ifndef NDEBUG
      // if (!check_simplex("multinomial_log_double", "Probabilites parameter", theta))
      //   return lp;
      // #endif
      for (unsigned int i = 0; i < ns.size(); ++i)
	lp += multiply_log(ns[i], theta[i]);
      return lp;
    }

    // double quadratic_form(vector_d& x, matrix_d& A) const {
    //   return (x.transpose() * (A.selfadjointView<Eigen::Lower>() * x))(0);
    // }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    template <bool propto, bool jacobian, bool lp_only>
    double log_prob_grad(vector<double>& params_r,
                         vector<double>& gradient,
                         ostream* pstream = 0) const {

      double lp = 0.0;
      paleon::Timer timer_total(1), timer_dss(K), timer_dgc(K), timer_exp(K), timer_inv(K), timer_drho(K);

      timer_total.tic(0);

      //
      // unpack model parameters and constrain
      //

      vector<double> vec_params_r;
      vector<int>    vec_params_i;
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r[i]);
      stan::io::reader<double> in(vec_params_r, vec_params_i);

      vector_d eta(K), eta_lja(K), eta_ja(K), eta_dj(K);
      for (int i=0; i<K; i++)
        eta[i] = lub_transform(in.scalar(), 0.1, 10.0, eta_lja[i], eta_ja[i], eta_dj[i]);

      if (jacobian)
        for (int i=0; i<K; i++)
          lp += eta_lja[i];

      vector_d rho(K), rho_lja(K), rho_ja(K), rho_dj(K);
      for (int i=0; i<K; i++)
        rho[i] = lub_transform(in.scalar(), 1.0e-6, 1.0, rho_lja[i], rho_ja[i], rho_dj[i]);

      if (jacobian)
        for (int i=0; i<K; i++)
          lp += rho_lja[i];

      vector_d mu = in.vector_constrain(K);

      vector<vector_d> alpha;
      for (int k = 0; k<K; ++k) {
        if (jacobian)
          alpha.push_back(in.vector_constrain(N_knots,lp));
        else
          alpha.push_back(in.vector_constrain(N_knots));
      }

      vector<vector_d> g;
      for (int k = 0; k<K; ++k) {
        if (jacobian)
          g.push_back(in.vector_constrain(N,lp));
        else
          g.push_back(in.vector_constrain(N));
      }

      //
      // compute log probability
      //

      double lp_thrd[K];

      vector<matrix_d> C(K);
      vector<matrix_d> c(K);
      vector<matrix_d> ct(K);
      vector<Eigen::LLT<matrix_d> > llt(K);

      matrix_d exp_g(N,K);
      matrix_d r(N,K);

      vector<matrix_d> C_inv(K);
      vector<vector_d> mu_g(K);
      vector<vector_d> sigma2(K);
      vector<vector_d> ci_Cinv_ci(K);
      vector<vector_d> cCalpha(K);

    
      lp += normal_log_double(mu, 0, 5);

      std::cout << "K = " << K <<  std::endl;

      #pragma omp parallel for
      for (int k=0; k<K; ++k) {
        double const eta2 = eta[k] * eta[k];
	double const rhoinv = 1.0 / rho[k];

	sigma2[k].resize(N);
	ci_Cinv_ci[k].resize(N);

	timer_exp.tic(k);
        C[k] = eta2 * (-rhoinv * d_knots.array()).exp().matrix();

        // XXX: 6%
        c[k] = eta2 * (-rhoinv * d_inter.array()).exp().matrix();
	timer_exp.toc(k);

        llt[k]  = C[k].llt();
	lp_thrd[k] = multi_normal_cholesky_log_double(alpha[k], zeros, llt[k].matrixL());

	timer_inv.tic(k);
        cCalpha[k] = c[k] * llt[k].solve(alpha[k]);
	C_inv[k] = C[k].inverse(); // XXX
	timer_inv.toc(k);

	mu_g[k] = mu[k] * ones + cCalpha[k] - Q * (QT * cCalpha[k]);

	for (int i=0; i<N; ++i) {
	  vector_d c_i = c[k].row(i);
          // XXX: 15%
	  ci_Cinv_ci[k][i] = (c_i.transpose() * C_inv[k] * c_i)(0);
	  sigma2[k][i] = eta2 + ci_Cinv_ci[k][i];

	  lp_thrd[k] += normal_log_double(g[k](i), mu_g[k][i], sqrt(sigma2[k][i]));
	}

	exp_g.col(k) = g[k].array().exp().matrix();
      }

     

      for (int k=0; k<K; ++k)
	lp += lp_thrd[k];
      
      vector_d sum_exp_g = exp_g.rowwise().sum();
      
      #pragma omp parallel for
      for (int k = 0; k < K; ++k){
        for (int i = 0; i < N; ++i){
          // r(i,k) = exp_g(i,k) / (1. + sum_exp_g(i));
	  r(i,k) = exp_g(i,k) / sum_exp_g(i);
	}
      }
      
      for (int i = 0; i < N; ++i) {
        //r(i,K) = 1. / (1. + sum_exp_g(i));
        if (r.row(i).sum() > 0.0)
          lp += multinomial_log_double(y[i], r.row(i));
      }
      
      if (lp_only)
        return lp;

      //
      // compute gradient
      //

      fill(gradient.begin(), gradient.end(), 0.0);

      #pragma omp parallel for
      for (int k=0; k<K; ++k) {

        double const rho2inv = 1.0 / (rho[k] * rho[k]);
        double const etainv  = 1.0 / eta[k];

	timer_drho.tic(k);
        vector_d const Cinv_alpha = llt[k].solve(alpha[k]);
        row_vector_d const alphaT = alpha[k].transpose();
        row_vector_d const alphaT_Cinv = Cinv_alpha.transpose();

        // XXX: 12%
        matrix_d const c_Cinv = llt[k].solve(c[k].transpose()).transpose();

        matrix_d const dCdrho = rho2inv * (C[k].array() * d_knots.array()).matrix();
        matrix_d const dcdrho = rho2inv * (c[k].array() * d_inter.array()).matrix();
	vector_d const dCdrho_Cinv_alpha = dCdrho * Cinv_alpha;
	timer_drho.toc(k);

	// mu
        gradient[2*K+k] -= mu[k] / 25;

        // eta
        gradient[k] += etainv * ( - N_knots + alphaT * Cinv_alpha ) * eta_ja[k] + eta_dj[k];

        // rho
        gradient[K+k] += 0.5 * (-(llt[k].solve(dCdrho)).trace()
                                 + alphaT_Cinv * dCdrho * Cinv_alpha) * rho_ja[k] + rho_dj[k];

        // mvn wrt alphas
        for (int i=0; i<N_knots; i++)
          gradient[3*K+N_knots*k+i] -= Cinv_alpha[i];

	timer_dgc.tic(k);
	row_vector_d const QT_c_Cinv = QT * c_Cinv;
	matrix_d const m_c_Cinv = (c_Cinv - Q * QT_c_Cinv).transpose();
	vector_d const dgcdrho = c_Cinv * (dCdrho * Cinv_alpha) - Q * (QT_c_Cinv * dCdrho_Cinv_alpha) - dcdrho * Cinv_alpha + Q * (QT * (dcdrho * Cinv_alpha));
	timer_dgc.toc(k);

	timer_dss.tic(k);
        // XXX: 4%
	vector_d const rs_cCinv_dcdrho = (c_Cinv.array() * dcdrho.array()).matrix().rowwise().sum();
        // XXX: 12%
	vector_d const rs_foo1 = ((dcdrho * C_inv[k]).array() * c[k].array()).matrix().rowwise().sum();
        // XXX: 12%
	vector_d const rs_foo2 = ((c[k] * (C_inv[k] * dCdrho * C_inv[k])).array() * c[k].array()).matrix().rowwise().sum();
	timer_dss.toc(k);

	// mu normal
	for (int i=0; i<N; i++) {

	  double const gc     = (g[k](i) - mu_g[k](i));
	  double const ss     = sigma2[k](i);
	  double const ssinv  = 1./ss;
	  double const gn     = gc * ssinv;
	  double const ss2    = ss * ss;
	  double const ss2inv = 1./ss2;

	  // XXX
	  double dssdrho = rs_foo1[i] - rs_foo2[i] + rs_cCinv_dcdrho[i];

	  // eta
	  gradient[k] += ((eta[k] + etainv * ci_Cinv_ci[k][i]) * gc*gc * ss2inv - etainv) * eta_ja[k];
	  gradient[2*K+k] += gn;

	  // g normal
	  gradient[3*K+K*N_knots+k*N+i] -= gn;

	  // rho normal
	  gradient[K+k] += -0.5 * gc * (2 * dgcdrho[i] * ss - gc * dssdrho) * ss2inv * rho_ja[k];
	  gradient[K+k] += -0.5 * ssinv * dssdrho * rho_ja[k];

	  for (int j=0; j<N_knots; j++)
	    gradient[3*K+k*N_knots+j] += gn * m_c_Cinv(j,i);

	}

	// partial of multinom
	for (int m=0; m<K; m++) {
	  for (int i=0; i<N; i++) {
	    //double const sumgp1 = 1. + sum_exp_g(i);
	    double const sumgp1 = sum_exp_g(i);
	    double const sumgp1inv2 = 1./(sumgp1*sumgp1);

	    double drdg;
            // if (m == K-1)
            //   drdg = -exp_g(i,k) * sumgp1inv2;
            // else if (m == k)
            //   drdg = exp_g(i,k) * (sumgp1 - exp_g(i,k)) * sumgp1inv2;
            // else if (m != k && m < K-1)
            //   drdg = -exp_g(i,k) * exp_g(i,m) * sumgp1inv2;
            if (m == k)
              drdg = exp_g(i,k) * (sumgp1 - exp_g(i,k)) * sumgp1inv2;
            else if (m != k)
              drdg = -exp_g(i,k) * exp_g(i,m) * sumgp1inv2;

	    gradient[3*K+K*N_knots+k*N+i] += y[i][m] / r(i,m) * drdg;
          }
        }
      }

      timer_total.toc(0);

      // cout << endl << scientific;
      // tinv.echo  ("    > inv          ");
      // tqf.echo   ("    > qf           ");
      // tda.echo   ("    > da           ");
      cout << scientific << endl;
      timer_exp.echo  (" > exp  ");
      timer_inv.echo  (" > inv  ");
      timer_drho.echo (" > drho ");
      timer_dgc.echo  (" > dgc  ");
      timer_dss.echo  (" > dss  ");
      timer_total.echo("=> lpg  ");

      return lp;
    }

    template <bool propto, bool jacobian>
    double log_prob(vector<double>& params_r,
                    vector<int>& params_i,
                    std::ostream* msgs = 0) const {

      vector<double> gradient;
      return log_prob_grad<propto, jacobian, true>(params_r, gradient, msgs);
    }

    void get_param_names(std::vector<std::string>& names__) const {
      names__.resize(0);
      names__.push_back("eta");
      names__.push_back("rho");
      names__.push_back("mu");
      names__.push_back("alpha");
      names__.push_back("g");
    }

    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
      dimss__.resize(0);
      std::vector<size_t> dims__;
      dims__.resize(0);
      dims__.push_back(K);
      dimss__.push_back(dims__);
      dims__.resize(0);
      dims__.push_back(K);
      dimss__.push_back(dims__);
      dims__.resize(0);
      dims__.push_back(K);
      dimss__.push_back(dims__);
      dims__.resize(0);
      dims__.push_back(K);
      dims__.push_back(N_knots);
      dimss__.push_back(dims__);
      dims__.resize(0);
      dims__.push_back(K);
      dims__.push_back(N);
      dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
      vars__.resize(0);
      stan::io::reader<double> in__(params_r__,params_i__);
      static const char* function__ = "veg_pp::write_array(%1%)";
      (void) function__; // dummy call to supress warning
      // read-transform, write parameters
      vector_d eta = in__.vector_lub_constrain(0.10000000000000001,sqrt(100),K);
      vector_d rho = in__.vector_lub_constrain(1.0e-6,1.0,K);
      vector_d mu = in__.vector_constrain(K);
      vector<vector_d> alpha;
      size_t dim_alpha_0__ = K;
      for (size_t k_0__ = 0; k_0__ < dim_alpha_0__; ++k_0__) {
	alpha.push_back(in__.vector_constrain(N_knots));
      }
      vector<vector_d> g;
      size_t dim_g_0__ = K;
      for (size_t k_0__ = 0; k_0__ < dim_g_0__; ++k_0__) {
	g.push_back(in__.vector_constrain(N));
      }
      for (int k_0__ = 0; k_0__ < K; ++k_0__) {
	vars__.push_back(eta[k_0__]);
      }
      for (int k_0__ = 0; k_0__ < K; ++k_0__) {
	vars__.push_back(rho[k_0__]);
      }
      for (int k_0__ = 0; k_0__ < K; ++k_0__) {
	vars__.push_back(mu[k_0__]);
      }
      for (int k_1__ = 0; k_1__ < N_knots; ++k_1__) {
	for (int k_0__ = 0; k_0__ < K; ++k_0__) {
	  vars__.push_back(alpha[k_0__][k_1__]);
	}
      }
      for (int k_1__ = 0; k_1__ < N; ++k_1__) {
	for (int k_0__ = 0; k_0__ < K; ++k_0__) {
	  vars__.push_back(g[k_0__][k_1__]);
	}
      }

      if (!include_tparams__) return;
      // declare and define transformed parameters
      double lp__ = 0.0;
      (void) lp__; // dummy call to supress warning
      stan::math::accumulator<double> lp_accum__;



      // validate transformed parameters

      // write transformed parameters

      if (!include_gqs__) return;
      // declare and define generated quantities


      // validate generated quantities

      // write generated quantities
    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }


    void write_csv_header(std::ostream& o__) const {
      stan::io::csv_writer writer__(o__);
      for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	writer__.comma();
	o__ << "eta" << '.' << k_0__;
      }
      for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	writer__.comma();
	o__ << "rho" << '.' << k_0__;
      }
      for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	writer__.comma();
	o__ << "mu" << '.' << k_0__;
      }
      for (int k_1__ = 1; k_1__ <= N_knots; ++k_1__) {
	for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	  writer__.comma();
	  o__ << "alpha" << '.' << k_0__ << '.' << k_1__;
	}
      }
      for (int k_1__ = 1; k_1__ <= N; ++k_1__) {
	for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	  writer__.comma();
	  o__ << "g" << '.' << k_0__ << '.' << k_1__;
	}
      }
      writer__.newline();
    }

    template <typename RNG>
    void write_csv(RNG& base_rng__,
                   std::vector<double>& params_r__,
                   std::vector<int>& params_i__,
                   std::ostream& o__,
                   std::ostream* pstream__ = 0) const {
      stan::io::reader<double> in__(params_r__,params_i__);
      stan::io::csv_writer writer__(o__);
      static const char* function__ = "veg_pp::write_csv(%1%)";
      (void) function__; // dummy call to supress warning
      // read-transform, write parameters
      vector_d eta = in__.vector_lub_constrain(0.10000000000000001,sqrt(100),K);
      writer__.write(eta);
      vector_d rho = in__.vector_lub_constrain(1.0e-6,1.0,K);
      writer__.write(rho);
      vector_d mu = in__.vector_constrain(K);
      writer__.write(mu);
      vector<vector_d> alpha;
      size_t dim_alpha_0__ = K;
      for (size_t k_0__ = 0; k_0__ < dim_alpha_0__; ++k_0__) {
	alpha.push_back(in__.vector_constrain(N_knots));
	writer__.write(alpha[k_0__]);
      }
      vector<vector_d> g;
      size_t dim_g_0__ = K;
      for (size_t k_0__ = 0; k_0__ < dim_g_0__; ++k_0__) {
	g.push_back(in__.vector_constrain(N));
	writer__.write(g[k_0__]);
      }

      // declare, define and validate transformed parameters
      double lp__ = 0.0;
      (void) lp__; // dummy call to supress warning
      stan::math::accumulator<double> lp_accum__;




      // write transformed parameters

      // declare and define generated quantities


      // validate generated quantities

      // write generated quantities
      writer__.newline();
    }

    template <typename RNG>
    void write_csv(RNG& base_rng,
                   Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                   std::ostream& o,
                   std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<int> params_i_vec;  // dummy
      write_csv(base_rng, params_r_vec, params_i_vec, o, pstream);
    }

    static std::string model_name() {
      return "veg_pp_model";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
      std::stringstream param_name_stream__;
      for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	param_name_stream__.str(std::string());
	param_name_stream__ << "eta" << '.' << k_0__;
	param_names__.push_back(param_name_stream__.str());
      }
      for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	param_name_stream__.str(std::string());
	param_name_stream__ << "rho" << '.' << k_0__;
	param_names__.push_back(param_name_stream__.str());
      }
      for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	param_name_stream__.str(std::string());
	param_name_stream__ << "mu" << '.' << k_0__;
	param_names__.push_back(param_name_stream__.str());
      }
      for (int k_1__ = 1; k_1__ <= N_knots; ++k_1__) {
	for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	  param_name_stream__.str(std::string());
	  param_name_stream__ << "alpha" << '.' << k_0__ << '.' << k_1__;
	  param_names__.push_back(param_name_stream__.str());
	}
      }
      for (int k_1__ = 1; k_1__ <= N; ++k_1__) {
	for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	  param_name_stream__.str(std::string());
	  param_name_stream__ << "g" << '.' << k_0__ << '.' << k_1__;
	  param_names__.push_back(param_name_stream__.str());
	}
      }

      if (!include_gqs__ && !include_tparams__) return;

      if (!include_gqs__) return;
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
      std::stringstream param_name_stream__;
      for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	param_name_stream__.str(std::string());
	param_name_stream__ << "eta" << '.' << k_0__;
	param_names__.push_back(param_name_stream__.str());
      }
      for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	param_name_stream__.str(std::string());
	param_name_stream__ << "rho" << '.' << k_0__;
	param_names__.push_back(param_name_stream__.str());
      }
      for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	param_name_stream__.str(std::string());
	param_name_stream__ << "mu" << '.' << k_0__;
	param_names__.push_back(param_name_stream__.str());
      }
      for (int k_1__ = 1; k_1__ <= N_knots; ++k_1__) {
	for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	  param_name_stream__.str(std::string());
	  param_name_stream__ << "alpha" << '.' << k_0__ << '.' << k_1__;
	  param_names__.push_back(param_name_stream__.str());
	}
      }
      for (int k_1__ = 1; k_1__ <= N; ++k_1__) {
	for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
	  param_name_stream__.str(std::string());
	  param_name_stream__ << "g" << '.' << k_0__ << '.' << k_1__;
	  param_names__.push_back(param_name_stream__.str());
	}
      }

      if (!include_gqs__ && !include_tparams__) return;

      if (!include_gqs__) return;
    }

  }; // model

} // namespace

typedef veg_pp::veg_pp_model stan_model;

namespace stan {
  namespace model {

    template <bool propto, bool jacobian_adjust_transform>
    double log_prob_grad(const stan_model& model,
                         std::vector<double>& params_r,
                         std::vector<int>& params_i,
                         std::vector<double>& gradient,
                         std::ostream* msgs = 0) {

      gradient.resize(params_r.size());
      return model.log_prob_grad<propto, jacobian_adjust_transform, false>(params_r, gradient, msgs);

    }

    void gradient(const stan_model& model,
                  const Eigen::Matrix<double, Eigen::Dynamic, 1>& x,
                  double& f,
                  Eigen::Matrix<double, Eigen::Dynamic, 1>& grad_f,
                  std::ostream* msgs = 0) {

      vector<double> params_r(x.rows());
      vector<double> grad(x.rows());

      for (size_t i = 0; i < params_r.size(); i++)
        params_r[i] = x[i];
      f = model.log_prob_grad<true, true, false>(params_r, grad, msgs);
      for (size_t i = 0; i < params_r.size(); i++)
        grad_f[i] = grad[i];
    }

  }
}
