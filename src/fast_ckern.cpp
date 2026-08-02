#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix build_vm_kernel_cpp(
    const NumericVector& data,
    const NumericVector& grid,
    const double kappa,
    const int n_threads) {
  const int n_data = data.size();
  const int n_grid = grid.size();
  NumericMatrix kernel(n_grid, n_data);

#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) schedule(static)
#endif
  for (int j = 0; j < n_data; ++j) {
    for (int g = 0; g < n_grid; ++g) {
      kernel(g, j) = std::exp(kappa * std::cos(grid[g] - data[j]));
    }
  }
  return kernel;
}

// [[Rcpp::export]]
NumericVector ckern_overlap_batch_cpp(
    const IntegerMatrix& samples,
    const NumericMatrix& kernel,
    const IntegerVector& left_index,
    const NumericVector& interpolation_weight,
    const int n_group1,
    const int n_threads) {
  const int n_reps = samples.nrow();
  const int n_total = samples.ncol();
  const int n_group2 = n_total - n_group1;
  const int n_grid = kernel.nrow();
  NumericVector result(n_reps);

#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) schedule(static)
#endif
  for (int r = 0; r < n_reps; ++r) {
    std::vector<double> density1(n_grid, 0.0);
    std::vector<double> density2(n_grid, 0.0);

    for (int i = 0; i < n_group1; ++i) {
      const int event = samples(r, i) - 1;
      for (int g = 0; g < n_grid; ++g) {
        density1[g] += kernel(g, event);
      }
    }
    for (int i = n_group1; i < n_total; ++i) {
      const int event = samples(r, i) - 1;
      for (int g = 0; g < n_grid; ++g) {
        density2[g] += kernel(g, event);
      }
    }

    const double scale1 = 1.0 / static_cast<double>(n_group1);
    const double scale2 = 1.0 / static_cast<double>(n_group2);
    for (int g = 0; g < n_grid; ++g) {
      density1[g] *= scale1;
      density2[g] *= scale2;
    }

    double overlap1 = 0.0;
    for (int i = 0; i < n_group1; ++i) {
      const int event = samples(r, i) - 1;
      const int left = left_index[event];
      const double alpha = interpolation_weight[event];
      const double f = density1[left] +
        alpha * (density1[left + 1] - density1[left]);
      const double g = density2[left] +
        alpha * (density2[left + 1] - density2[left]);
      overlap1 += std::min(1.0, g / f);
    }

    double overlap2 = 0.0;
    for (int i = n_group1; i < n_total; ++i) {
      const int event = samples(r, i) - 1;
      const int left = left_index[event];
      const double alpha = interpolation_weight[event];
      const double f = density1[left] +
        alpha * (density1[left + 1] - density1[left]);
      const double g = density2[left] +
        alpha * (density2[left + 1] - density2[left]);
      overlap2 += std::min(1.0, f / g);
    }

    result[r] = 0.5 * (
      overlap1 / static_cast<double>(n_group1) +
      overlap2 / static_cast<double>(n_group2)
    );
  }
  return result;
}
