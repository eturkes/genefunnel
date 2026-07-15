// Assisted-by: OpenAI Codex.

/*
 *    This file is part of GeneFunnel.
 *    Copyright (C) 2025-2026  Emir Turkes, UK DRI at UCL
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
 *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
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

double finalize_score(
  const std::size_t observed_size,
  const long double sum_values,
  const std::size_t below_mean_count,
  const long double below_mean_sum,
  const long double value_scale,
  const int gene_set,
  const int column
) {
  if (!std::isfinite(sum_values) ||
      !std::isfinite(below_mean_sum) ||
      !std::isfinite(value_scale) ||
      value_scale <= 0.0L ||
      below_mean_count >= observed_size) {
    stop_non_finite_result(gene_set, column);
  }

  // Since total absolute deviation equals twice the deviation below the
  // mean, the normative formula is equivalently:
  // ((n - 1 - l) * sum(x) + n * sum(x[x < mean(x)])) / (n - 1),
  // where l is the number below the mean. For non-negative inputs every
  // term is non-negative, avoiding catastrophic subtraction near zero.
  const long double denominator = observed_size - 1;
  const std::size_t total_weight =
    observed_size - 1 - below_mean_count;
  long double total_term = total_weight * sum_values;
  long double below_mean_term = observed_size * below_mean_sum;
  if (std::isfinite(total_term) && std::isfinite(below_mean_term)) {
    total_term /= denominator;
    below_mean_term /= denominator;
  } else {
    // Form bounded coefficients when multiplying either sum by n would
    // overflow despite a potentially representable final score.
    total_term = (total_weight / denominator) * sum_values;
    below_mean_term =
      (observed_size / denominator) * below_mean_sum;
  }
  const long double score_in_units = total_term + below_mean_term;
  const long double score_extended = score_in_units * value_scale;
  const long double scale = std::max(
    std::fabs(total_term),
    std::fabs(below_mean_term)
  ) * value_scale;
  if (!std::isfinite(score_extended) || !std::isfinite(scale)) {
    stop_non_finite_result(gene_set, column);
  }

  const long double tolerance = score_tolerance(observed_size, scale);
  if (score_extended < -tolerance) {
    stop(
      "Internal GeneFunnel arithmetic produced a materially negative "
      "score for gene set %d, column %d.",
      gene_set,
      column
    );
  }

  const double score = static_cast<double>(score_extended);
  if (!std::isfinite(score)) {
    stop_non_finite_result(gene_set, column);
  }

  // Clamp negative roundoff only. Small positive scores may be genuine.
  return score <= 0.0 ? 0.0 : score;
}

double score_with_scaled_sum(
  const std::size_t observed_size,
  const std::vector<double>& positive_values,
  const int gene_set,
  const int column
) {
  if (positive_values.empty()) {
    stop_non_finite_result(gene_set, column);
  }

  const double maximum = *std::max_element(
    positive_values.begin(),
    positive_values.end()
  );
  long double normalized_sum = 0.0L;
  for (const double value : positive_values) {
    normalized_sum += static_cast<long double>(value) / maximum;
  }

  const long double normalized_mean = normalized_sum / observed_size;
  std::size_t below_mean_count =
    observed_size - positive_values.size();
  long double below_mean_sum = 0.0L;
  for (const double value : positive_values) {
    const long double normalized_value =
      static_cast<long double>(value) / maximum;
    if (normalized_value < normalized_mean) {
      ++below_mean_count;
      below_mean_sum += normalized_value;
    }
  }

  return finalize_score(
    observed_size,
    normalized_sum,
    below_mean_count,
    below_mean_sum,
    maximum,
    gene_set,
    column
  );
}

template <typename ValueAt>
NumericMatrix calculate_dense_scores(
  const int n_rows,
  const int n_columns,
  const List& gene_indices,
  const ValueAt& value_at
) {
  const int n_gene_sets = gene_indices.size();

  NumericMatrix scores(n_gene_sets, n_columns);

  for (int j = 0; j < n_columns; ++j) {
    for (int i = 0; i < n_gene_sets; ++i) {
      IntegerVector indices = gene_indices[i];

      if (indices.size() < 2) {
        scores(i, j) = NA_REAL;
        continue;
      }

      std::vector<double> observed;
      observed.reserve(indices.size());

      for (R_xlen_t k = 0; k < indices.size(); ++k) {
        int index = indices[k];
        if (index == NA_INTEGER || index < 1 || index > n_rows) {
          stop("Internal gene-set index is invalid.");
        }

        double value = value_at(index - 1, j);
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
        if (value < 0.0) {
          stop(
            "Internal GeneFunnel arithmetic produced a materially negative "
            "score for gene set %d, column %d.",
            i + 1,
            j + 1
          );
        }

        observed.push_back(value);
      }

      const std::size_t observed_size = observed.size();
      if (observed_size < 2) {
        scores(i, j) = NA_REAL;
        continue;
      }

      long double sum_values = 0.0L;
      for (const double value : observed) {
        sum_values += static_cast<long double>(value);
      }
      if (!std::isfinite(sum_values)) {
        std::vector<double> positive_values;
        positive_values.reserve(observed_size);
        for (const double value : observed) {
          if (value > 0.0) {
            positive_values.push_back(value);
          }
        }
        scores(i, j) = score_with_scaled_sum(
          observed_size,
          positive_values,
          i + 1,
          j + 1
        );
        continue;
      }

      const long double mean_values = sum_values / observed_size;
      // On platforms where long double equals double, a positive subnormal
      // mean can round to zero. Exact zeros remain below that mathematical
      // mean; positive observed values do not.
      const bool mean_underflowed =
        sum_values > 0.0L && mean_values == 0.0L;
      std::size_t below_mean_count = 0;
      long double below_mean_sum = 0.0L;
      for (const double value : observed) {
        const bool below_mean = mean_underflowed
          ? value == 0.0
          : static_cast<long double>(value) < mean_values;
        if (below_mean) {
          ++below_mean_count;
          below_mean_sum += static_cast<long double>(value);
        }
      }

      scores(i, j) = finalize_score(
        observed_size,
        sum_values,
        below_mean_count,
        below_mean_sum,
        1.0L,
        i + 1,
        j + 1
      );
    }
  }

  return scores;
}

struct SparseMembership {
  int gene_set;
  R_xlen_t position;
};

struct SparseEntry {
  R_xlen_t position;
  double value;
};

NumericMatrix calculate_sparse_scores(
  const arma::sp_mat& orig_mat,
  const List& gene_indices
) {
  const int n_rows = static_cast<int>(orig_mat.n_rows);
  const int n_columns = static_cast<int>(orig_mat.n_cols);
  const int n_gene_sets = gene_indices.size();
  NumericMatrix scores(n_gene_sets, n_columns);

  // Invert memberships once per bounded native chunk. Each sparse column can
  // then be streamed once; implicit zeros require neither lookup nor storage.
  std::vector<std::vector<SparseMembership>> memberships_by_row(n_rows);
  std::vector<std::size_t> gene_set_sizes(n_gene_sets);
  for (int i = 0; i < n_gene_sets; ++i) {
    IntegerVector indices = gene_indices[i];
    gene_set_sizes[i] = indices.size();
    for (R_xlen_t k = 0; k < indices.size(); ++k) {
      const int index = indices[k];
      if (index == NA_INTEGER || index < 1 || index > n_rows) {
        stop("Internal gene-set index is invalid.");
      }
      memberships_by_row[index - 1].push_back({i, k});
    }
  }

  std::vector<std::vector<SparseEntry>> entries_by_gene_set(n_gene_sets);
  for (int j = 0; j < n_columns; ++j) {
    for (auto& entries : entries_by_gene_set) {
      entries.clear();
    }

    for (arma::sp_mat::const_col_iterator entry = orig_mat.begin_col(j);
         entry != orig_mat.end_col(j);
         ++entry) {
      const int row = static_cast<int>(entry.row());
      const double value = *entry;
      const auto& memberships = memberships_by_row[row];
      if (memberships.empty()) {
        continue;
      }
      if (std::isinf(value)) {
        stop(
          "Internal input contains an infinite value at matrix row %d, "
          "column %d, gene set %d.",
          row + 1,
          j + 1,
          memberships.front().gene_set + 1
        );
      }

      for (const SparseMembership& membership : memberships) {
        entries_by_gene_set[membership.gene_set].push_back({
          membership.position,
          value
        });
      }
    }

    for (int i = 0; i < n_gene_sets; ++i) {
      const std::size_t gene_set_size = gene_set_sizes[i];
      if (gene_set_size < 2) {
        scores(i, j) = NA_REAL;
        continue;
      }

      auto& entries = entries_by_gene_set[i];
      std::sort(
        entries.begin(),
        entries.end(),
        [](const SparseEntry& left, const SparseEntry& right) {
          return left.position < right.position;
        }
      );

      std::size_t observed_size = gene_set_size;
      std::size_t positive_count = 0;
      long double sum_values = 0.0L;
      for (const SparseEntry& entry : entries) {
        if (std::isnan(entry.value)) {
          --observed_size;
          continue;
        }
        if (entry.value < 0.0) {
          stop(
            "Internal GeneFunnel arithmetic produced a materially negative "
            "score for gene set %d, column %d.",
            i + 1,
            j + 1
          );
        }
        if (entry.value > 0.0) {
          ++positive_count;
          sum_values += static_cast<long double>(entry.value);
        }
      }

      if (observed_size < 2) {
        scores(i, j) = NA_REAL;
        continue;
      }

      if (!std::isfinite(sum_values)) {
        std::vector<double> positive_values;
        positive_values.reserve(positive_count);
        for (const SparseEntry& entry : entries) {
          if (!std::isnan(entry.value) && entry.value > 0.0) {
            positive_values.push_back(entry.value);
          }
        }
        scores(i, j) = score_with_scaled_sum(
          observed_size,
          positive_values,
          i + 1,
          j + 1
        );
        continue;
      }

      const long double mean_values = sum_values / observed_size;
      // A positive mathematical mean can underflow when long double has the
      // same range as double. Implicit and stored zeros are still below it.
      const bool positive_mean = sum_values > 0.0L;
      std::size_t below_mean_count = positive_mean
        ? observed_size - positive_count
        : 0;
      long double below_mean_sum = 0.0L;
      for (const SparseEntry& entry : entries) {
        if (!std::isnan(entry.value) &&
            entry.value > 0.0 &&
            static_cast<long double>(entry.value) < mean_values) {
          ++below_mean_count;
          below_mean_sum += static_cast<long double>(entry.value);
        }
      }

      scores(i, j) = finalize_score(
        observed_size,
        sum_values,
        below_mean_count,
        below_mean_sum,
        1.0L,
        i + 1,
        j + 1
      );
    }
  }

  return scores;
}

}  // namespace

// [[Rcpp::export]]
NumericMatrix calculateScoresDense(
  const NumericMatrix& orig_mat,
  const List& gene_indices
) {
  return calculate_dense_scores(
    orig_mat.nrow(),
    orig_mat.ncol(),
    gene_indices,
    [&orig_mat](const int row, const int column) {
      return orig_mat(row, column);
    }
  );
}

// [[Rcpp::export]]
NumericMatrix calculateScoresSparse(
  const arma::sp_mat& orig_mat,
  const List& gene_indices
) {
  return calculate_sparse_scores(orig_mat, gene_indices);
}
