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

#include <dolphin/common/common_base.h>

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matexpr/mat_arith.h>
#include <light_mat/matexpr/mat_emath.h>
#include <light_mat/matexpr/mat_pred.h>
#include <light_mat/matexpr/mat_cast.h>

namespace dolphin
{
	using lmat::dimension;
	using lmat::matrix_shape;

	using lmat::IMatrixXpr;
	using lmat::IEWiseMatrix;
	using lmat::IRegularMatrix;

	using lmat::copy;
	using lmat::fill;
	using lmat::zero;
	using lmat::copy_from;

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

	using lmat::in_;
	using lmat::out_;
	using lmat::in_out_;

	using lmat::ewise;
	using lmat::percol;

	// import math functions

	using lmat::max;
	using lmat::min;
	using lmat::cond;
	using lmat::clamp;
	using lmat::fma;

	using lmat::abs;
	using lmat::sqr;
	using lmat::cube;
	using lmat::cbrt;
	using lmat::rcp;
	using lmat::rsqrt;

	using lmat::floor;
	using lmat::ceil;
	using lmat::round;
	using lmat::trunc;

	using lmat::signbit;
	using lmat::isfinite;
	using lmat::isinf;
	using lmat::isnan;

	using lmat::exp;
	using lmat::log;
	using lmat::log10;
	using lmat::xlogx;
	using lmat::xlogy;

	using lmat::exp2;
	using lmat::log2;
	using lmat::expm1;
	using lmat::log1p;

	using lmat::sin;
	using lmat::cos;
	using lmat::tan;

	using lmat::asin;
	using lmat::acos;
	using lmat::atan;
	using lmat::atan2;

	using lmat::sinh;
	using lmat::cosh;
	using lmat::tanh;

	using lmat::asinh;
	using lmat::acosh;
	using lmat::atanh;

}

#endif
