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
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

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

      vec idx_values(indices.size());
      for (R_xlen_t k = 0; k < indices.size(); ++k) {
        int index = indices[k];
        if (index == NA_INTEGER || index < 1 ||
            static_cast<uword>(index) > orig_mat.n_rows) {
          stop("Internal gene-set index is invalid.");
        }
        idx_values[k] = orig_mat(index - 1, j);
      }

      double sum_values = sum(idx_values);
      double var_values = sum(abs(idx_values - mean(idx_values)));

      uword size = idx_values.size();
      double factor = static_cast<double>(size) / (2.0 * (size - 1));
      double score = sum_values - (var_values * factor);

      double epsilon = 1e-9;
      if (fabs(score) < epsilon) {
        score = 0.0;
      }

      mat(i, j) = score;
    }
  }

  return mat;
}
