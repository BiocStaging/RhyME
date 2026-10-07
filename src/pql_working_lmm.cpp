// [[Rcpp::depends(RcppEigen, BH)]]
#include "pql_working_lmm_numeric.h"
#include <boost/math/tools/minima.hpp>

namespace {

// R supplies columns in the order: intercept, covariates, sine, cosine.
Design read_design(const NumericMatrix& design) {
    const int n = design.nrow(), p = design.ncol();
    Design X;
    X.reserve(p);
    for (int j=0; j<p; ++j)
        X.emplace_back(design.begin()+j*n,design.begin()+(j+1)*n);
    return X;
}

struct WorkingFit {
    bool valid = false;
    double objective = std::numeric_limits<double>::infinity();
    double log_ratio = NA_REAL, d = NA_REAL, tau2 = NA_REAL, quadratic = NA_REAL;
    Coefficients beta;
    Covariance covariance;
    explicit WorkingFit(int p) : beta(p,0), covariance(p*p,0) {}
};

double score_statistic(const Eigen::Vector2d& score, Eigen::Matrix2d information) {
    information = (0.5*(information+information.transpose())).eval();
    Eigen::LLT<Eigen::Matrix2d> factor(information);
    const double scale = std::max(information.diagonal().cwiseAbs().maxCoeff(),1e-12);
    if (!information.allFinite() || factor.info()!=Eigen::Success ||
        !(factor.matrixL().toDenseMatrix().diagonal().array().square()>1e-12*scale).all())
        stop("Singular efficient score information.");
    return std::max(0.0,score.dot(factor.solve(score)));
}

struct WorkingProblem {
    const Design& X;
    const std::vector<double>& z;
    const std::vector<double>& w0;
    const IntegerVector& id;
    const std::vector<double>& group_w0;
    const int p;
    bool fixed_only;
    bool use_prior;
    double prior_mean, prior_sd;
    double logdet_A = 0;
    int evaluations = 0;
    const double d_lower = 1.0 / (1.0 + std::exp(15.0));
    const double d_upper = 1.0 / (1.0 + std::exp(-5.0));

    WorkingProblem(const Design& X_,
                   const std::vector<double>& z_, const std::vector<double>& w0_,
                   const IntegerVector& id_, const std::vector<double>& gw_,
                   bool null_model_ = false, bool use_prior_ = false,
                   double prior_mean_ = 0, double prior_sd_ = 1, bool fixed_only_ = false)
        : X(X_), z(z_), w0(w0_), id(id_), group_w0(gw_),
          p(X_.size()-(null_model_ ? 2 : 0)),
          fixed_only(fixed_only_), use_prior(use_prior_), prior_mean(prior_mean_), prior_sd(prior_sd_) {
        for (double w : w0) logdet_A -= std::log(w);
    }

    double prior_penalty(double d) const {
        const double standardized = (std::log1p(-d)-std::log(d)-prior_mean)/prior_sd;
        return 0.5*standardized*standardized;
    }

    double profile_scale(double quadratic, double n, double lo, double hi) const {
        const double ml_d = clamp_value(quadratic/n,lo,hi);
        if (!use_prior || hi-lo<1e-15) return ml_d;
        // With a prior, optimize log(phi); the penalty needs no Jacobian.
        const double lower = std::log1p(-hi)-std::log(hi);
        const double upper = std::log1p(-lo)-std::log(lo);
        auto objective = [&](double u) {
            const double d = 1/(1+std::exp(u));
            return 0.5*(n*std::log(d)+quadratic/d)+prior_penalty(d);
        };
        std::vector<double> grid, values;
        for (int i=0;i<=16;++i) grid.push_back(lower+(upper-lower)*i/16);
        grid.push_back(clamp_value(prior_mean,lower,upper));
        grid.push_back(std::log1p(-ml_d)-std::log(ml_d));
        std::sort(grid.begin(),grid.end());
        grid.erase(std::unique(grid.begin(),grid.end()),grid.end());
        double best_u = grid[0], best = objective(best_u);
        for (double u : grid) {
            values.push_back(objective(u));
            if (values.back()<best) {best=values.back(); best_u=u;}
        }
        for (std::size_t i=1;i+1<grid.size();++i) {
            if (!(values[i]<=values[i-1] && values[i]<=values[i+1] &&
                  (values[i]<values[i-1] || values[i]<values[i+1]))) continue;
            boost::math::uintmax_t iterations=100;
            const auto minimum=boost::math::tools::brent_find_minima(
                objective,grid[i-1],grid[i+1],26,iterations);
            if (iterations>=100 || !std::isfinite(minimum.second))
                stop("Prior scale Brent refinement failed.");
            if (minimum.second<best) {best=minimum.second; best_u=minimum.first;}
        }
        return clamp_value(1/(1+std::exp(best_u)),lo,hi);
    }

    WorkingFit evaluate(double log_ratio) {
        ++evaluations;
        WorkingFit fit(p);
        if (fixed_only) log_ratio = -std::numeric_limits<double>::infinity();
        fit.log_ratio = log_ratio;
        const double ratio = std::exp(log_ratio); // tau^2 / d
        const double lo = fixed_only ? d_lower : std::max(d_lower, kTau2Lower / ratio);
        const double hi = fixed_only ? d_upper : std::min(d_upper, kTau2Upper / ratio);
        if (lo > hi * (1.0 + 1e-12)) return fit;
        Covariance inverse_information(p*p,0);
        if (!gls_fit(log_ratio,X,p,z,w0,id,group_w0,fit.beta,inverse_information)) return fit;

        double logdet_H = logdet_A, quadratic = 0;
        std::vector<double> residual(z.size()), group_mean(group_w0.size(), 0.0);
        for (std::size_t i = 0; i < z.size(); ++i) {
            residual[i] = residual_at(z[i],X,fit.beta,i);
            if (fixed_only) quadratic += w0[i]*residual[i]*residual[i];
            else group_mean[id[i]] += w0[i]*residual[i];
        }
        if (!fixed_only) {
            // Center within subjects to avoid cancellation in the residual quadratic.
            for (std::size_t j = 0; j < group_w0.size(); ++j) {
                if (group_w0[j] > 0) group_mean[j] /= group_w0[j];
                quadratic += group_w0[j] * group_mean[j] * group_mean[j] /
                             (1.0 + ratio * group_w0[j]);
                logdet_H += std::log1p(ratio * group_w0[j]);
            }
            for (std::size_t i = 0; i < z.size(); ++i) {
                const double centered = residual[i] - group_mean[id[i]];
                quadratic += w0[i] * centered * centered;
            }
        }
        const double n = z.size();
        fit.d = profile_scale(quadratic, n, lo, std::max(lo, hi));
        fit.tau2 = ratio * fit.d;
        fit.quadratic = quadratic;
        // log|V| = n log(d) + log|H|.
        fit.objective = 0.5 * (n * std::log(fit.d) + logdet_H + quadratic / fit.d);
        if (use_prior) fit.objective += prior_penalty(fit.d);
        for (int k = 0; k < p*p; ++k) fit.covariance[k] = fit.d * inverse_information[k];
        fit.valid = std::isfinite(fit.objective);
        return fit;
    }

    WorkingFit optimize() {
        if (fixed_only) return evaluate(-std::numeric_limits<double>::infinity());
        // Bounds for tau^2 / d implied by the precision and variance limits.
        const double lower = std::log(kTau2Lower / d_upper);
        const double upper = std::log(kTau2Upper / d_lower);
        const int intervals = 40;
        std::vector<double> grid(intervals + 1), values(intervals + 1);
        WorkingFit best(p);
        for (int j = 0; j <= intervals; ++j) {
            grid[j] = lower + (upper - lower) * j / intervals;
            WorkingFit fit = evaluate(grid[j]);
            values[j] = fit.objective;
            if (fit.valid && fit.objective < best.objective) best = fit;
        }
        // Refine each local grid minimum; the objective may be multimodal.
        for (int j = 1; j < intervals; ++j) {
            if (!(values[j] <= values[j-1] && values[j] <= values[j+1] &&
                  (values[j] < values[j-1] || values[j] < values[j+1]))) continue;
            boost::math::uintmax_t iterations = 100;
            const auto minimum = boost::math::tools::brent_find_minima(
                [&](double x) { return evaluate(x).objective; },
                grid[j-1], grid[j+1], 26, iterations);
            const WorkingFit fit = evaluate(minimum.first);
            if (iterations >= 100 || !fit.valid)
                stop("Working-LMM Brent refinement failed.");
            if (fit.objective < best.objective) best = fit;
        }
        return best;
    }
};

List working_fit_list(const WorkingFit& fit, int evaluations) {
    const int columns = fit.beta.size();
    NumericMatrix covariance(columns,columns);
    for (int j=0; j<columns; ++j) for (int k=0; k<columns; ++k)
        covariance(j,k) = fit.valid ? fit.covariance[columns*j+k] : NA_REAL;
    return List::create(_["valid"] = fit.valid, _["objective"] = fit.objective,
        _["d"] = fit.d, _["tau2"] = fit.tau2, _["log_ratio"] = fit.log_ratio,
        _["beta"] = wrap(fit.beta), _["covariance"] = covariance,
        _["quadratic_H"] = fit.quadratic, _["evaluations"] = evaluations);
}
} // namespace

// Fit a frozen working model for numerical checks.
// w0 = 1 / diag(A), not the precision-scaled final working weights.
// [[Rcpp::export]]
List pql_working_lmm_step_cpp(NumericVector z, NumericMatrix design,
        IntegerVector subject, NumericVector w0,
        double log_ratio = NA_REAL, bool null_model = false,
        bool use_prior = false, double prior_mean = 0, double prior_sd = 1,
        bool fixed_only = false) {
    const int n = z.size();
    if (design.nrow()!=n || subject.size()!=n || w0.size()!=n)
        stop("Invalid working-model dimensions.");
    if (use_prior && (!std::isfinite(prior_mean) || !std::isfinite(prior_sd) || prior_sd<=0))
        stop("Prior requires finite mean and positive finite SD.");
    int groups = 0;
    const Design X = read_design(design);
    for (double value : design)
        if (!std::isfinite(value)) stop("Nonfinite working-model design.");
    for (int i=0; i<n; ++i) {
        if (!std::isfinite(z[i]) ||
            !std::isfinite(w0[i]) || w0[i]<=0 || subject[i]<0 || subject[i]>=n)
            stop("Invalid working-model input.");
        groups = std::max(groups, subject[i]+1);
    }
    std::vector<double> weights = as<std::vector<double>>(w0), response = as<std::vector<double>>(z);
    std::vector<double> gw(groups, 0);
    for (int i=0; i<n; ++i) gw[subject[i]] += w0[i];
    WorkingProblem problem(X, response, weights, subject, gw, null_model,
                           use_prior, prior_mean, prior_sd, fixed_only);
    WorkingFit fit = std::isfinite(log_ratio) ? problem.evaluate(log_ratio) : problem.optimize();
    return working_fit_list(fit, problem.evaluations);
}

// Detailed working quantities are returned only for numerical checks.
// Response, design and prior are validated by the R fitting entry point.
// [[Rcpp::export]]
List pql_working_lmm_cpp(NumericVector y, NumericMatrix design, IntegerVector subject,
        int max_iterations = 500, double tolerance = 1e-5,
        double epsilon = 1e-4, double tau_absolute_tolerance = 1e-7,
        double tau_relative_tolerance = 1e-4, double update_relaxation = 0.75,
        double working_response_limit = 1000.0, bool null_model = false,
        bool use_prior = false, double prior_mean = 0, double prior_sd = 1,
        bool fixed_only = false, bool diagnostics = false) {
    const int n = y.size();
    if (design.nrow()!=n || subject.size()!=n ||
        max_iterations<1 || tolerance<=0 || tau_absolute_tolerance<=0 ||
        tau_relative_tolerance<=0 || !(epsilon>0 && epsilon<0.5) ||
        !(update_relaxation>0 && update_relaxation<=1) ||
        !std::isfinite(working_response_limit) || working_response_limit<0)
        stop("Invalid PQL fit arguments.");
    y = clone(y);
    int groups = 0;
    double mean_y = 0;
    const Design X = read_design(design);
    for (int i=0; i<n; ++i) {
        if (subject[i]<0 || subject[i]>=n) stop("Invalid subject index.");
        y[i] = clamp_value(y[i], epsilon, 1-epsilon);
        mean_y += y[i]; groups = std::max(groups, subject[i]+1);
    }
    mean_y /= n;
    double variance = 0;
    for (double v : y) variance += (v-mean_y)*(v-mean_y);
    variance /= n-1;
    double log_phi = std::log(clamp_value(mean_y*(1-mean_y)/std::max(variance,1e-8)-1,
                                         std::exp(-5.0), std::exp(15.0)));
    double tau2 = fixed_only ? 0.0 : 0.01;
    const int nuisance = X.size()-2;
    const int p = null_model ? nuisance : X.size();
    Coefficients beta(p,0);
    beta[0] = std::log(mean_y/(1-mean_y));
    std::vector<double> u(groups,0), eta(n), mu(n), z(n), w0(n), gw(groups,0);
    int iterations = 0, evaluations = 0, capped = 0;
    bool converged = false;
    auto make_working = [&]() {
        std::fill(gw.begin(),gw.end(),0.0);
        for (int i=0; i<n; ++i) {
            eta[i] = beta[0]*X[0][i];
            for (int j=1; j<p; ++j) eta[i] += beta[j]*X[j][i];
            eta[i] += u[subject[i]];
            mu[i] = clamp_value(inv_logit(eta[i]),1e-7,1-1e-7);
            w0[i] = mu[i]*(1-mu[i]);
            z[i] = eta[i]+(y[i]-mu[i])/w0[i];
            gw[subject[i]] += w0[i];
        }
        capped += cap_extreme_working_response(z,working_response_limit);
    };
    WorkingFit step(p);
    for (int it=0; it<max_iterations; ++it) {
        iterations = it+1;
        make_working();
        WorkingProblem problem(X,z,w0,subject,gw,null_model,use_prior,prior_mean,prior_sd,fixed_only);
        step = problem.optimize();
        evaluations += problem.evaluations;
        if (!step.valid) stop("Invalid joint working-LMM variance update.");
        double next_log_phi = std::log1p(-step.d)-std::log(step.d);
        double next_tau2 = step.tau2;
        Coefficients next_beta = step.beta;
        std::vector<double> next_u(groups,0);
        if (!fixed_only) {
            for (int i=0; i<n; ++i) next_u[subject[i]] += w0[i] *
                residual_at(z[i],X,next_beta,i);
            const double ratio = next_tau2/step.d;
            for (int j=0; j<groups; ++j) next_u[j] *= ratio/(1+ratio*gw[j]);
            next_tau2 = damp_variance(next_tau2, tau2, update_relaxation);
            damp_updates(next_u, u, update_relaxation);
        }
        damp_updates(next_beta, beta, update_relaxation);
        next_log_phi = damp_update(next_log_phi, log_phi, update_relaxation);
        const bool stable =
            stable_parameter(next_log_phi, log_phi, tolerance) &&
            stable_variance(next_tau2, tau2, tau_absolute_tolerance, tau_relative_tolerance) &&
            stable_parameters(next_beta, beta, tolerance) &&
            stable_parameters(next_u, u, tolerance);
        beta=next_beta; log_phi=next_log_phi; tau2=next_tau2; u.swap(next_u);
        if (stable) {converged=true; break;}
    }
    make_working();
    const double d = 1/(1+std::exp(log_phi));
    std::vector<double> weight(n), group_weight(groups,0);
    for (int i=0; i<n; ++i) {weight[i]=w0[i]/d; group_weight[subject[i]]+=weight[i];}
    Covariance covariance(p*p,0);
    if (!gls_fit(std::log(tau2),X,p,z,weight,subject,group_weight,beta,covariance))
        stop("Singular final fixed-effect information.");
    double statistic = NA_REAL;
    if (null_model) {
        std::vector<double> residual(n);
        for (int i=0; i<n; ++i) residual[i] = residual_at(z[i],X,beta,i);
        const auto vi_residual = apply_v_inverse(residual,weight,subject,group_weight,tau2);
        Eigen::Vector2d score;
        Eigen::Matrix2d information;
        Eigen::MatrixXd cross(nuisance,2), inverse(nuisance,nuisance);
        for (int j=0; j<nuisance; ++j) for (int k=0; k<nuisance; ++k)
            inverse(j,k) = covariance[p*j+k];
        for (int j=0; j<2; ++j) {
            score[j] = dot_product(X[nuisance+j],vi_residual);
            const auto vi_x = apply_v_inverse(X[nuisance+j],weight,subject,group_weight,tau2);
            for (int k=0; k<nuisance; ++k) cross(k,j) = dot_product(X[k],vi_x);
            for (int k=0; k<2; ++k) information(k,j) = dot_product(X[nuisance+k],vi_x);
        }
        information -= cross.transpose()*inverse*cross;
        statistic = score_statistic(score,information);
    } else if (diagnostics) {
        const int s = nuisance, c = nuisance+1;
        const double ss = covariance[p*s+s], cc = covariance[p*c+c];
        const double sc = covariance[p*s+c];
        const double determinant = ss*cc-sc*sc;
        if (!(determinant>0)) stop("Singular sine/cosine covariance.");
        statistic = std::max(0.0,(cc*beta[s]*beta[s]-2*sc*beta[s]*beta[c]+
            ss*beta[c]*beta[c])/determinant);
    }

    if ((null_model || diagnostics) && !std::isfinite(statistic))
        stop("Nonfinite test statistic.");
    // Classify the final undamped proposal using the convergence tolerances.
    const double proposed_log_phi = std::log1p(-step.d)-std::log(step.d);
    const bool tau_upper = !fixed_only && stable_variance(step.tau2, kTau2Upper,
        tau_absolute_tolerance, tau_relative_tolerance);
    const bool tau_lower = !fixed_only && tau2<=kTau2Lower*(1+1e-8);
    const bool phi_bound = stable_parameter(proposed_log_phi, -5.0, tolerance) ||
        stable_parameter(proposed_log_phi, 15.0, tolerance);
    const char* status = tau_upper ? "boundary_upper" : phi_bound ? "precision_boundary" :
        converged ? (tau_lower ? "boundary_zero" : "converged") : "max_iterations";
    if (!diagnostics) {
        if (null_model)
            return List::create(_["log_phi"]=log_phi,_["score_chisq"]=statistic,
                _["p_value"]=R::pchisq(statistic,2,false,false),
                _["converged"]=converged && !tau_upper && !phi_bound,_["status"]=status);
        return List::create(_["b0"]=beta[0],_["b_sin"]=beta[nuisance],
            _["b_cos"]=beta[nuisance+1],_["log_phi"]=log_phi,
            _["tau"]=tau_lower ? 0.0 : std::sqrt(tau2),
            _["converged"]=converged && !tau_upper && !phi_bound,_["status"]=status);
    }
    NumericMatrix cov(p,p);
    for (int j=0; j<p; ++j) for (int k=0; k<p; ++k)
        cov(j,k) = covariance[p*j+k];
    List answer = List::create(_["b0"]=beta[0],_["b_sin"]=null_model ? NA_REAL : beta[nuisance],
        _["b_cos"]=null_model ? NA_REAL : beta[nuisance+1],
        _["log_phi"]=log_phi,_["tau"]=tau_lower ? 0.0 : std::sqrt(tau2),
        _["tau2_internal"]=tau2,_["d"]=d,_["iterations"]=iterations,
        _["converged"]=converged && !tau_upper && !phi_bound,_["status"]=status,
        _["wald_chisq"]=statistic,_["p_value"]=R::pchisq(statistic,2,false,false),
        _["u"]=u,_["covariance"]=cov,_["working_z"]=z,_["working_w0"]=w0,
        _["profile_evaluations"]=evaluations,_["capped_working_values"]=capped,
        _["beta"]=beta);
    if (null_model) {
        answer["score_chisq"] = statistic;
        answer["wald_chisq"] = R_NilValue;
    }
    return answer;
}
