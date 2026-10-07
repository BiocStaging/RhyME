#pragma once
// Numerical helpers for the random-intercept cosinor model.
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

using namespace Rcpp;

namespace {

constexpr double kTau2Lower = 1e-8;
constexpr double kTau2Upper = 10.0;

inline double clamp_value(double x, double lower, double upper) {
    return std::max(lower, std::min(upper, x));
}

inline double inv_logit(double x) {
    if (x >= 0.0) {
        const double e = std::exp(-x);
        return 1.0 / (1.0 + e);
    }
    const double e = std::exp(x);
    return e / (1.0 + e);
}

std::vector<double> apply_v_inverse(const std::vector<double>& value,
                                    const std::vector<double>& weight,
                                    const IntegerVector& subject,
                                    const std::vector<double>& group_weight,
                                    double tau2) {
    if (tau2==0) {
        std::vector<double> answer(value.size());
        for (std::size_t i=0; i<value.size(); ++i) answer[i] = weight[i]*value[i];
        return answer;
    }
    const int n_subjects = group_weight.size();
    std::vector<double> group_weighted_value(n_subjects, 0.0);
    for (R_xlen_t i = 0; i < subject.size(); ++i) {
        group_weighted_value[subject[i]] += weight[i] * value[i];
    }

    std::vector<double> answer(value.size());
    for (R_xlen_t i = 0; i < subject.size(); ++i) {
        const int group = subject[i];
        const double denominator = 1.0 + tau2 * group_weight[group];
        answer[i] = weight[i] * value[i] -
            tau2 * weight[i] * group_weighted_value[group] / denominator;
    }
    return answer;
}

double dot_product(const std::vector<double>& x,
                   const std::vector<double>& y) {
    return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}

int cap_extreme_working_response(std::vector<double>& value,
                                double limit) {
    if (!(limit > 0.0)) return 0;
    if (std::none_of(value.begin(), value.end(),
                     [limit](double x) { return std::abs(x) > limit; })) return 0;
    std::vector<double> sorted = value;
    const std::size_t middle = sorted.size() / 2;
    std::nth_element(
        sorted.begin(), sorted.begin() + middle, sorted.end()
    );
    double median = sorted[middle];
    if (sorted.size() % 2 == 0) {
        std::nth_element(
            sorted.begin(), sorted.begin() + middle - 1, sorted.end()
        );
        median = 0.5 * (median + sorted[middle - 1]);
    }
    int capped = 0;
    for (double& x : value)
        if (std::abs(x) > limit) { x = median; ++capped; }
    return capped;
}


inline bool stable_parameter(double next, double current, double tolerance) {
    return std::abs(next - current) <=
        tolerance * (1.0 + std::max(std::abs(next), std::abs(current)));
}

template<class Container>
bool stable_parameters(const Container& next, const Container& current,
                       double tolerance) {
    for (std::size_t j = 0; j < next.size(); ++j)
        if (!stable_parameter(next[j], current[j], tolerance)) return false;
    return true;
}

// The absolute tolerance handles variance estimates near zero.
inline bool stable_variance(double next, double current,
                            double absolute_tolerance, double relative_tolerance) {
    return std::abs(next - current) <=
        absolute_tolerance + relative_tolerance * std::max(next, current);
}

inline double damp_update(double next, double current, double relaxation) {
    return relaxation < 1.0 ? current + relaxation * (next - current) : next;
}

inline double damp_variance(double next, double current, double relaxation) {
    return relaxation < 1.0 ?
        std::exp(std::log(current) + relaxation * (std::log(next) - std::log(current))) :
        next;
}

template<class Container>
void damp_updates(Container& next, const Container& current, double relaxation) {
    if (relaxation < 1.0)
        for (std::size_t j = 0; j < next.size(); ++j)
            next[j] = current[j] + relaxation * (next[j] - current[j]);
}

using Design = std::vector<std::vector<double>>;
using Coefficients = std::vector<double>;
using Covariance = std::vector<double>; // row-major

bool solve_information(const Covariance& information, const Coefficients& rhs,
                       int p, Coefficients& beta, Covariance& covariance) {
    using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    const Eigen::Map<const RowMatrix> value(information.data(), p, p);
    if (!value.allFinite()) return false;
    Eigen::LLT<Eigen::MatrixXd> factor(value);
    if (factor.info() != Eigen::Success) return false;
    const double scale = std::max(value.diagonal().cwiseAbs().maxCoeff(), 1e-12);
    if (!(factor.matrixL().toDenseMatrix().diagonal().array().square() > 1e-12*scale).all())
        return false;
    const Eigen::VectorXd fitted = factor.solve(Eigen::Map<const Eigen::VectorXd>(rhs.data(), p));
    const Eigen::MatrixXd inverse = factor.solve(Eigen::MatrixXd::Identity(p, p));
    if (!fitted.allFinite() || !inverse.allFinite()) return false;
    std::copy(fitted.data(), fitted.data()+p, beta.begin());
    for (int row=0; row<p; ++row) for (int column=0; column<p; ++column)
        covariance[p*row+column] = 0.5*(inverse(row,column)+inverse(column,row));
    return true;
}

bool gls_fit(double log_tau2, const Design& design, int p,
             const std::vector<double>& response, const std::vector<double>& weight,
             const IntegerVector& subject, const std::vector<double>& group_weight,
             Coefficients& beta, Covariance& covariance) {
    const double tau2 = std::exp(log_tau2);
    Design vi_design(p);
    for (int column=0; column<p; ++column)
        vi_design[column] = apply_v_inverse(design[column],weight,subject,group_weight,tau2);
    const auto vi_response = apply_v_inverse(response,weight,subject,group_weight,tau2);
    Covariance information(p*p,0);
    Coefficients rhs(p,0);
    for (int row=0; row<p; ++row) {
        rhs[row] = dot_product(design[row],vi_response);
        for (int column=0; column<p; ++column)
            information[p*row+column] = dot_product(design[row],vi_design[column]);
    }
    return solve_information(information,rhs,p,beta,covariance);
}

double residual_at(double response, const Design& design,
                   const Coefficients& beta, int i) {
    double fitted = 0;
    for (std::size_t j=0; j<beta.size(); ++j) fitted += beta[j]*design[j][i];
    return response-fitted;
}
} // namespace
