/*
 *    This file is part of GeneFunnel.
 *    Copyright (C) 2025  Emir Turkes, UK DRI at UCL
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *    Emir Turkes can be contacted at emir.turkes@eturkes.com
 */

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

long double score_tolerance(
  const std::size_t observed_size,
  const long double scale
) {
  const long double operation_bound = 8.0L * observed_size + 16.0L;
  return operation_bound * std::numeric_limits<double>::epsilon() * scale;
}

void stop_non_finite_result(const int gene_set, const int column) {
  stop(
    "Internal GeneFunnel arithmetic produced a non-finite result for "
    "gene set %d, column %d.",
    gene_set,
    column
  );
}

}  // namespace

// [[Rcpp::export]]
NumericMatrix calculateScores(
  const arma::sp_mat& orig_mat, List gene_indices
) {
  int ncol_mat = orig_mat.n_cols;
  int nrow_list = gene_indices.size();

  NumericMatrix mat(nrow_list, ncol_mat);

  for (int j = 0; j < ncol_mat; ++j) {
    for (int i = 0; i < nrow_list; ++i) {
      IntegerVector indices = gene_indices[i];

      if (indices.size() < 2) {
        mat(i, j) = NA_REAL;
        continue;
      }

      std::vector<double> observed;
      observed.reserve(indices.size());

      for (R_xlen_t k = 0; k < indices.size(); ++k) {
        int index = indices[k];
        if (index == NA_INTEGER || index < 1 ||
            static_cast<uword>(index) > orig_mat.n_rows) {
          stop("Internal gene-set index is invalid.");
        }

        double value = orig_mat(index - 1, j);
        if (std::isnan(value)) {
          continue;
        }
        if (std::isinf(value)) {
          stop(
            "Internal input contains an infinite value at matrix row %d, "
            "column %d, gene set %d.",
            index,
            j + 1,
            i + 1
          );
        }

        observed.push_back(value);
      }

      const std::size_t observed_size = observed.size();
      if (observed_size < 2) {
        mat(i, j) = NA_REAL;
        continue;
      }

      long double sum_values = 0.0L;
      for (const double value : observed) {
        sum_values += static_cast<long double>(value);
      }
      if (!std::isfinite(sum_values)) {
        stop_non_finite_result(i + 1, j + 1);
      }

      const long double mean_values = sum_values / observed_size;
      std::size_t below_mean_count = 0;
      long double below_mean_sum = 0.0L;
      for (const double value : observed) {
        if (static_cast<long double>(value) < mean_values) {
          ++below_mean_count;
          below_mean_sum += static_cast<long double>(value);
        }
      }
      if (!std::isfinite(below_mean_sum) ||
          below_mean_count >= observed_size) {
        stop_non_finite_result(i + 1, j + 1);
      }

      // Since total absolute deviation equals twice the deviation below the
      // mean, the normative formula is equivalently:
      // ((n - 1 - l) * sum(x) + n * sum(x[x < mean(x)])) / (n - 1),
      // where l is the number below the mean. For non-negative inputs every
      // term is non-negative, avoiding catastrophic subtraction near zero.
      const long double total_term =
        (observed_size - 1 - below_mean_count) * sum_values;
      const long double below_mean_term =
        observed_size * below_mean_sum;
      const long double score_extended =
        (total_term + below_mean_term) / (observed_size - 1);
      const long double scale = std::max(
        std::fabs(total_term),
        std::fabs(below_mean_term)
      ) / (observed_size - 1);
      if (!std::isfinite(score_extended) || !std::isfinite(scale)) {
        stop_non_finite_result(i + 1, j + 1);
      }

      const long double tolerance = score_tolerance(
        observed_size,
        scale
      );
      if (score_extended < -tolerance) {
        stop(
          "Internal GeneFunnel arithmetic produced a materially negative "
          "score for gene set %d, column %d.",
          i + 1,
          j + 1
        );
      }

      const double score = static_cast<double>(score_extended);
      if (!std::isfinite(score)) {
        stop_non_finite_result(i + 1, j + 1);
      }

      // Clamp negative roundoff only. Small positive scores may be genuine.
      mat(i, j) = score <= 0.0 ? 0.0 : score;
    }
  }

  return mat;
}
