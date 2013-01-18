/**
 * @file import_lmat.h
 *
 * @brief Import useful matrix classes from Light-Matrix
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef DOLPHIN_IMPORT_LMAT_H_
#define DOLPHIN_IMPORT_LMAT_H_

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matexpr/mat_arith.h>
#include <light_mat/matexpr/mat_emath.h>
#include <light_mat/matexpr/mat_pred.h>
#include <light_mat/matexpr/mat_cast.h>

namespace dolphin
{
	using lmat::IMatrixXpr;
	using lmat::IEWiseMatrix;
	using lmat::IRegularMatrix;

	using lmat::dense_matrix;
	using lmat::dense_col;
	using lmat::dense_row;

	using lmat::ref_matrix;
	using lmat::ref_col;
	using lmat::ref_grid;

	using lmat::ref_block;
	using lmat::ref_grid;
	using lmat::step_col;
	using lmat::step_row;

	using lmat::meta::is_contiguous;
	using lmat::meta::is_percol_contiguous;
	using lmat::meta::supports_linear_index;

	using lmat::is_vector;
	using lmat::is_column;
	using lmat::is_row;
	using lmat::is_empty;
	using lmat::is_scalar;
	using lmat::is_square;

}

#endif
