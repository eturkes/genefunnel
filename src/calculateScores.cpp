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

double score_observed_values(
  const std::vector<double>& observed,
  const int gene_set,
  const int column
) {
  const std::size_t observed_size = observed.size();
  if (observed_size < 2) {
    return NA_REAL;
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
    return score_with_scaled_sum(
      observed_size,
      positive_values,
      gene_set,
      column
    );
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

  return finalize_score(
    observed_size,
    sum_values,
    below_mean_count,
    below_mean_sum,
    1.0L,
    gene_set,
    column
  );
}

// Error-free transforms provide a platform-independent double-double
// significand. A separate exponent keeps intermediate diagnostics normalized
// across the complete finite binary64 range.
struct SumError {
  double sum;
  double error;
};

SumError two_sum(const double left, const double right) {
  const double sum = left + right;
  const double right_virtual = sum - left;
  const double error =
    (left - (sum - right_virtual)) + (right - right_virtual);
  return {sum, error};
}

SumError two_product(const double left, const double right) {
  const double product = left * right;
  constexpr double splitter = 134217729.0;
  const double left_split = splitter * left;
  const double right_split = splitter * right;
  const double left_high = left_split - (left_split - left);
  const double right_high = right_split - (right_split - right);
  const double left_low = left - left_high;
  const double right_low = right - right_high;
  const double error =
    ((left_high * right_high - product) +
      left_high * right_low + left_low * right_high) +
    left_low * right_low;
  return {product, error};
}

struct ScaledValue {
  double high;
  double low;
  int exponent;
};

double component_error_allowance(const std::size_t observed_size) {
  const double unit_roundoff =
    std::numeric_limits<double>::epsilon() / 2.0;
  return 64.0 * (static_cast<double>(observed_size) + 1.0) *
    unit_roundoff;
}

ScaledValue normalize_scaled(
  double high,
  double low = 0.0,
  int exponent = 0
) {
  SumError combined = two_sum(high, low);
  high = combined.sum;
  low = combined.error;
  if (high == 0.0) {
    if (low == 0.0) {
      return {0.0, 0.0, 0};
    }
    high = low;
    low = 0.0;
  }
  if (!std::isfinite(high) || !std::isfinite(low)) {
    stop("Internal component arithmetic could not be certified.");
  }

  int adjustment = 0;
  std::frexp(std::fabs(high), &adjustment);
  high = std::scalbn(high, -adjustment);
  low = std::scalbn(low, -adjustment);
  exponent += adjustment;

  combined = two_sum(high, low);
  high = combined.sum;
  low = combined.error;
  if (std::fabs(high) >= 1.0) {
    high *= 0.5;
    low *= 0.5;
    ++exponent;
  } else if (std::fabs(high) < 0.5) {
    high *= 2.0;
    low *= 2.0;
    --exponent;
  }
  return {high, low, exponent};
}

ScaledValue scaled_from_double(const double value) {
  if (!std::isfinite(value)) {
    stop("Internal component arithmetic received a non-finite value.");
  }
  return normalize_scaled(value);
}

ScaledValue negate_scaled(const ScaledValue& value) {
  return {-value.high, -value.low, value.exponent};
}

ScaledValue add_scaled(
  const ScaledValue& left,
  const ScaledValue& right
) {
  if (left.high == 0.0) {
    return right;
  }
  if (right.high == 0.0) {
    return left;
  }

  const int exponent = std::max(left.exponent, right.exponent);
  const int left_shift = left.exponent - exponent;
  const int right_shift = right.exponent - exponent;
  const double left_high = std::scalbn(left.high, left_shift);
  const double left_low = std::scalbn(left.low, left_shift);
  const double right_high = std::scalbn(right.high, right_shift);
  const double right_low = std::scalbn(right.low, right_shift);

  const SumError leading = two_sum(left_high, right_high);
  const SumError trailing = two_sum(left_low, right_low);
  const SumError middle = two_sum(leading.error, trailing.sum);
  const SumError result = two_sum(leading.sum, middle.sum);
  const double low = result.error + middle.error + trailing.error;
  return normalize_scaled(result.sum, low, exponent);
}

int compare_scaled(
  const ScaledValue& left,
  const ScaledValue& right
) {
  if (left.high < 0.0 || right.high < 0.0) {
    stop("Internal component comparison requires non-negative values.");
  }
  if (left.high == 0.0 && right.high == 0.0) {
    return 0;
  }
  if (left.high == 0.0) {
    return -1;
  }
  if (right.high == 0.0) {
    return 1;
  }
  if (left.exponent != right.exponent) {
    return left.exponent < right.exponent ? -1 : 1;
  }
  if (left.high != right.high) {
    return left.high < right.high ? -1 : 1;
  }
  if (left.low != right.low) {
    return left.low < right.low ? -1 : 1;
  }
  return 0;
}

ScaledValue subtract_scaled(
  const ScaledValue& left,
  const ScaledValue& right
) {
  if (compare_scaled(left, right) < 0) {
    stop("Internal component subtraction would be negative.");
  }
  return add_scaled(left, negate_scaled(right));
}

ScaledValue multiply_scaled(
  const ScaledValue& left,
  const ScaledValue& right
) {
  if (left.high == 0.0 || right.high == 0.0) {
    return scaled_from_double(0.0);
  }
  const SumError leading = two_product(left.high, right.high);
  const double cross =
    left.high * right.low + left.low * right.high;
  const double trailing =
    leading.error + cross + left.low * right.low;
  const SumError combined = two_sum(leading.sum, trailing);
  return normalize_scaled(
    combined.sum,
    combined.error,
    left.exponent + right.exponent
  );
}

ScaledValue multiply_scaled_double(
  const ScaledValue& value,
  const double multiplier
) {
  if (!std::isfinite(multiplier) || multiplier < 0.0) {
    stop("Internal component multiplier is invalid.");
  }
  return multiply_scaled(value, scaled_from_double(multiplier));
}

ScaledValue divide_scaled(
  const ScaledValue& numerator,
  const ScaledValue& denominator
) {
  if (numerator.high < 0.0 || denominator.high <= 0.0) {
    stop("Internal component division is invalid.");
  }
  if (numerator.high == 0.0) {
    return scaled_from_double(0.0);
  }

  ScaledValue quotient = normalize_scaled(
    numerator.high / denominator.high,
    0.0,
    numerator.exponent - denominator.exponent
  );
  for (int iteration = 0; iteration < 3; ++iteration) {
    const ScaledValue product = multiply_scaled(denominator, quotient);
    const int comparison = compare_scaled(numerator, product);
    if (comparison == 0) {
      break;
    }
    const ScaledValue residual = comparison > 0
      ? subtract_scaled(numerator, product)
      : subtract_scaled(product, numerator);
    const ScaledValue correction = normalize_scaled(
      residual.high / denominator.high,
      0.0,
      residual.exponent - denominator.exponent
    );
    quotient = comparison > 0
      ? add_scaled(quotient, correction)
      : subtract_scaled(quotient, correction);
  }
  return quotient;
}

ScaledValue divide_scaled_double(
  const ScaledValue& value,
  const double divisor
) {
  if (!std::isfinite(divisor) || divisor <= 0.0) {
    stop("Internal component divisor is invalid.");
  }
  return divide_scaled(value, scaled_from_double(divisor));
}

ScaledValue absolute_scaled_difference(
  const ScaledValue& left,
  const ScaledValue& right
) {
  return compare_scaled(left, right) >= 0
    ? subtract_scaled(left, right)
    : subtract_scaled(right, left);
}

struct ScaledPair {
  double mantissa;
  int exponent;
};

ScaledPair scaled_pair(const ScaledValue& value) {
  ScaledValue normalized = normalize_scaled(
    value.high,
    value.low,
    value.exponent
  );
  if (normalized.high == 0.0) {
    return {0.0, 0};
  }
  double mantissa = normalized.high + normalized.low;
  if (std::fabs(mantissa) >= 1.0) {
    mantissa *= 0.5;
    ++normalized.exponent;
  } else if (std::fabs(mantissa) < 0.5) {
    mantissa *= 2.0;
    --normalized.exponent;
  }
  return {mantissa, normalized.exponent};
}

double scaled_as_double(const ScaledValue& value) {
  const ScaledPair pair = scaled_pair(value);
  const double ordinary = std::scalbn(pair.mantissa, pair.exponent);
  if (pair.mantissa != 0.0 &&
      (!std::isfinite(ordinary) || ordinary == 0.0)) {
    return NA_REAL;
  }
  return ordinary;
}

ScaledValue sum_scaled_values(const std::vector<double>& values) {
  ScaledValue result = scaled_from_double(0.0);
  for (const double value : values) {
    result = add_scaled(result, scaled_from_double(value));
  }
  return result;
}

// Fixed native mirror of the public result schema. Diagnostic matrices exist
// only in calculateComponents*(); calculateScores*() retains its original
// score-only allocation path.
struct ComponentMatrices {
  NumericMatrix score;
  NumericMatrix observed_sum;
  NumericMatrix penalty;
  NumericMatrix balance;
  IntegerMatrix effective_size;
  NumericMatrix observed_fraction;
  CharacterMatrix semantic_status;
  CharacterMatrix observed_sum_status;
  CharacterMatrix penalty_status;
  CharacterMatrix balance_status;
  CharacterMatrix conditioning_status;
  NumericMatrix observed_sum_mantissa;
  IntegerMatrix observed_sum_exponent;
  NumericMatrix penalty_mantissa;
  IntegerMatrix penalty_exponent;
  NumericMatrix balance_mantissa;
  IntegerMatrix balance_exponent;

  ComponentMatrices(const int n_gene_sets, const int n_columns) :
    score(n_gene_sets, n_columns),
    observed_sum(n_gene_sets, n_columns),
    penalty(n_gene_sets, n_columns),
    balance(n_gene_sets, n_columns),
    effective_size(n_gene_sets, n_columns),
    observed_fraction(n_gene_sets, n_columns),
    semantic_status(n_gene_sets, n_columns),
    observed_sum_status(n_gene_sets, n_columns),
    penalty_status(n_gene_sets, n_columns),
    balance_status(n_gene_sets, n_columns),
    conditioning_status(n_gene_sets, n_columns),
    observed_sum_mantissa(n_gene_sets, n_columns),
    observed_sum_exponent(n_gene_sets, n_columns),
    penalty_mantissa(n_gene_sets, n_columns),
    penalty_exponent(n_gene_sets, n_columns),
    balance_mantissa(n_gene_sets, n_columns),
    balance_exponent(n_gene_sets, n_columns) {
    std::fill(score.begin(), score.end(), NA_REAL);
    std::fill(observed_sum.begin(), observed_sum.end(), NA_REAL);
    std::fill(penalty.begin(), penalty.end(), NA_REAL);
    std::fill(balance.begin(), balance.end(), NA_REAL);
    std::fill(observed_fraction.begin(), observed_fraction.end(), NA_REAL);
    std::fill(
      observed_sum_mantissa.begin(),
      observed_sum_mantissa.end(),
      NA_REAL
    );
    std::fill(
      observed_sum_exponent.begin(),
      observed_sum_exponent.end(),
      NA_INTEGER
    );
    std::fill(
      penalty_mantissa.begin(),
      penalty_mantissa.end(),
      NA_REAL
    );
    std::fill(
      penalty_exponent.begin(),
      penalty_exponent.end(),
      NA_INTEGER
    );
    std::fill(
      balance_mantissa.begin(),
      balance_mantissa.end(),
      NA_REAL
    );
    std::fill(
      balance_exponent.begin(),
      balance_exponent.end(),
      NA_INTEGER
    );
  }

  List as_list() const {
    return List::create(
      _["score"] = score,
      _["observed_sum"] = observed_sum,
      _["penalty"] = penalty,
      _["balance"] = balance,
      _["effective_size"] = effective_size,
      _["observed_fraction"] = observed_fraction,
      _["status"] = List::create(
        _["semantic"] = semantic_status,
        _["observed_sum"] = observed_sum_status,
        _["penalty"] = penalty_status,
        _["balance"] = balance_status,
        _["conditioning"] = conditioning_status
      ),
      _["scaled"] = List::create(
        _["observed_sum"] = List::create(
          _["mantissa"] = observed_sum_mantissa,
          _["exponent"] = observed_sum_exponent
        ),
        _["penalty"] = List::create(
          _["mantissa"] = penalty_mantissa,
          _["exponent"] = penalty_exponent
        ),
        _["balance"] = List::create(
          _["mantissa"] = balance_mantissa,
          _["exponent"] = balance_exponent
        )
      )
    );
  }
};

bool store_component_value(
  const ScaledValue& value,
  NumericMatrix& ordinary,
  CharacterMatrix& status,
  NumericMatrix& mantissa,
  IntegerMatrix& exponent,
  const int gene_set,
  const int column
) {
  if (value.high < 0.0) {
    stop("Internal component arithmetic produced a negative diagnostic.");
  }
  const double ordinary_value = scaled_as_double(value);
  if (!NumericVector::is_na(ordinary_value)) {
    ordinary(gene_set, column) = ordinary_value;
    status(gene_set, column) = "ordinary";
    return true;
  }

  const ScaledPair pair = scaled_pair(value);
  ordinary(gene_set, column) = NA_REAL;
  status(gene_set, column) = "scaled";
  mantissa(gene_set, column) = pair.mantissa;
  exponent(gene_set, column) = pair.exponent;
  return false;
}

void store_unavailable_component(
  NumericMatrix& ordinary,
  CharacterMatrix& status,
  const int gene_set,
  const int column
) {
  ordinary(gene_set, column) = NA_REAL;
  status(gene_set, column) = "unavailable";
}

void store_component_cell(
  ComponentMatrices& result,
  const int gene_set,
  const int column,
  const std::size_t matched_size,
  const std::size_t observed_size,
  const std::vector<double>& positive_values,
  const double score
) {
  if (matched_size == 0 || observed_size > matched_size ||
      positive_values.size() > observed_size) {
    stop("Internal component support counts are invalid.");
  }
  result.score(gene_set, column) = score;
  result.effective_size(gene_set, column) =
    static_cast<int>(observed_size);
  result.observed_fraction(gene_set, column) =
    static_cast<double>(observed_size) / matched_size;

  const ScaledValue observed_sum = sum_scaled_values(positive_values);
  const bool sum_is_ordinary = store_component_value(
    observed_sum,
    result.observed_sum,
    result.observed_sum_status,
    result.observed_sum_mantissa,
    result.observed_sum_exponent,
    gene_set,
    column
  );

  if (observed_size < 2) {
    result.semantic_status(gene_set, column) = "too_few_observed";
    store_unavailable_component(
      result.penalty,
      result.penalty_status,
      gene_set,
      column
    );
    store_unavailable_component(
      result.balance,
      result.balance_status,
      gene_set,
      column
    );
    result.conditioning_status(gene_set, column) = "not_applicable";
    return;
  }

  if (observed_sum.high == 0.0) {
    result.semantic_status(gene_set, column) = "zero_total";
    store_component_value(
      scaled_from_double(0.0),
      result.penalty,
      result.penalty_status,
      result.penalty_mantissa,
      result.penalty_exponent,
      gene_set,
      column
    );
    store_unavailable_component(
      result.balance,
      result.balance_status,
      gene_set,
      column
    );
    result.conditioning_status(gene_set, column) = "not_applicable";
    return;
  }

  result.semantic_status(gene_set, column) = "scoreable";
  const ScaledValue center = divide_scaled_double(
    observed_sum,
    static_cast<double>(observed_size)
  );
  const std::size_t zero_count = observed_size - positive_values.size();
  ScaledValue deviation_sum = multiply_scaled_double(
    center,
    static_cast<double>(zero_count)
  );
  std::size_t below_mean_count = zero_count;
  ScaledValue below_mean_sum = scaled_from_double(0.0);
  for (const double value : positive_values) {
    const ScaledValue scaled_value = scaled_from_double(value);
    deviation_sum = add_scaled(
      deviation_sum,
      absolute_scaled_difference(scaled_value, center)
    );
    if (compare_scaled(scaled_value, center) < 0) {
      ++below_mean_count;
      below_mean_sum = add_scaled(below_mean_sum, scaled_value);
    }
  }
  if (below_mean_count >= observed_size) {
    stop("Internal component mean classification is invalid.");
  }

  const ScaledValue penalty = divide_scaled_double(
    multiply_scaled_double(
      deviation_sum,
      static_cast<double>(observed_size)
    ),
    2.0 * static_cast<double>(observed_size - 1)
  );
  const ScaledValue scaled_score = divide_scaled_double(
    add_scaled(
      multiply_scaled_double(
        observed_sum,
        static_cast<double>(
          observed_size - 1 - below_mean_count
        )
      ),
      multiply_scaled_double(
        below_mean_sum,
        static_cast<double>(observed_size)
      )
    ),
    static_cast<double>(observed_size - 1)
  );
  ScaledValue balance = divide_scaled(scaled_score, observed_sum);
  const ScaledValue one = scaled_from_double(1.0);
  if (compare_scaled(balance, one) > 0) {
    const ScaledValue excess = subtract_scaled(balance, one);
    const ScaledValue rounding_budget = scaled_from_double(
      8.0 * component_error_allowance(observed_size)
    );
    if (compare_scaled(excess, rounding_budget) > 0) {
      stop("Internal component balance exceeded its mathematical bound.");
    }
    balance = one;
  }

  const bool penalty_is_ordinary = store_component_value(
    penalty,
    result.penalty,
    result.penalty_status,
    result.penalty_mantissa,
    result.penalty_exponent,
    gene_set,
    column
  );
  const bool balance_is_ordinary = store_component_value(
    balance,
    result.balance,
    result.balance_status,
    result.balance_mantissa,
    result.balance_exponent,
    gene_set,
    column
  );

  bool safe = sum_is_ordinary && penalty_is_ordinary &&
    balance_is_ordinary && std::isfinite(score) && score > 0.0 &&
    score >= std::numeric_limits<double>::min() &&
    scaled_score.high > 0.0;
  if (safe) {
    const ScaledValue kappa = divide_scaled(
      add_scaled(observed_sum, penalty),
      scaled_score
    );
    const double ordinary_kappa = scaled_as_double(kappa);
    const double allowance = component_error_allowance(observed_size);
    safe = std::isfinite(ordinary_kappa) && allowance < 1.0 &&
      allowance * ordinary_kappa <= std::ldexp(1.0, -20);
  }
  result.conditioning_status(gene_set, column) = safe
    ? "safe"
    : "ill_conditioned";
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

      scores(i, j) = score_observed_values(observed, i + 1, j + 1);
    }
  }

  return scores;
}

List calculate_dense_components(
  const NumericMatrix& orig_mat,
  const List& gene_indices
) {
  const int n_rows = orig_mat.nrow();
  const int n_columns = orig_mat.ncol();
  const int n_gene_sets = gene_indices.size();
  ComponentMatrices result(n_gene_sets, n_columns);

  // Read each requested matrix cell once; subsequent numerical passes use the
  // bounded observed-value buffer, never the input matrix again.
  for (int j = 0; j < n_columns; ++j) {
    for (int i = 0; i < n_gene_sets; ++i) {
      IntegerVector indices = gene_indices[i];
      std::vector<double> observed;
      observed.reserve(indices.size());

      for (R_xlen_t k = 0; k < indices.size(); ++k) {
        const int index = indices[k];
        if (index == NA_INTEGER || index < 1 || index > n_rows) {
          stop("Internal gene-set index is invalid.");
        }

        const double value = orig_mat(index - 1, j);
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

      const double score = score_observed_values(
        observed,
        i + 1,
        j + 1
      );
      std::vector<double> positive_values;
      positive_values.reserve(observed.size());
      for (const double value : observed) {
        if (value > 0.0) {
          positive_values.push_back(value);
        }
      }
      store_component_cell(
        result,
        i,
        j,
        indices.size(),
        observed.size(),
        positive_values,
        score
      );
    }
  }

  return result.as_list();
}

struct SparseMembership {
  int gene_set;
  R_xlen_t position;
};

struct SparseEntry {
  R_xlen_t position;
  double value;
};

double score_sparse_entries(
  const std::size_t gene_set_size,
  const std::vector<SparseEntry>& entries,
  const int gene_set,
  const int column
) {
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
        gene_set,
        column
      );
    }
    if (entry.value > 0.0) {
      ++positive_count;
      sum_values += static_cast<long double>(entry.value);
    }
  }

  if (observed_size < 2) {
    return NA_REAL;
  }

  if (!std::isfinite(sum_values)) {
    std::vector<double> positive_values;
    positive_values.reserve(positive_count);
    for (const SparseEntry& entry : entries) {
      if (!std::isnan(entry.value) && entry.value > 0.0) {
        positive_values.push_back(entry.value);
      }
    }
    return score_with_scaled_sum(
      observed_size,
      positive_values,
      gene_set,
      column
    );
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

  return finalize_score(
    observed_size,
    sum_values,
    below_mean_count,
    below_mean_sum,
    1.0L,
    gene_set,
    column
  );
}

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

      scores(i, j) = score_sparse_entries(
        gene_set_size,
        entries,
        i + 1,
        j + 1
      );
    }
  }

  return scores;
}

List calculate_sparse_components(
  const arma::sp_mat& orig_mat,
  const List& gene_indices
) {
  const int n_rows = static_cast<int>(orig_mat.n_rows);
  const int n_columns = static_cast<int>(orig_mat.n_cols);
  const int n_gene_sets = gene_indices.size();
  ComponentMatrices result(n_gene_sets, n_columns);

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

    // Stream each sparse column once; implicit zeros remain represented by
    // observed_size - positive_values.size().
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
      auto& entries = entries_by_gene_set[i];
      std::sort(
        entries.begin(),
        entries.end(),
        [](const SparseEntry& left, const SparseEntry& right) {
          return left.position < right.position;
        }
      );

      const std::size_t gene_set_size = gene_set_sizes[i];
      const double score = score_sparse_entries(
        gene_set_size,
        entries,
        i + 1,
        j + 1
      );
      std::size_t observed_size = gene_set_size;
      std::vector<double> positive_values;
      positive_values.reserve(entries.size());
      for (const SparseEntry& entry : entries) {
        if (std::isnan(entry.value)) {
          --observed_size;
        } else if (entry.value > 0.0) {
          positive_values.push_back(entry.value);
        }
      }
      store_component_cell(
        result,
        i,
        j,
        gene_set_size,
        observed_size,
        positive_values,
        score
      );
    }
  }

  return result.as_list();
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

// [[Rcpp::export]]
List calculateComponentsDense(
  const NumericMatrix& orig_mat,
  const List& gene_indices
) {
  return calculate_dense_components(orig_mat, gene_indices);
}

// [[Rcpp::export]]
List calculateComponentsSparse(
  const arma::sp_mat& orig_mat,
  const List& gene_indices
) {
  return calculate_sparse_components(orig_mat, gene_indices);
}
