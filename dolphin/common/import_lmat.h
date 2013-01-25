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
#include <light_mat/matrix/matrix_asvec.h>
#include <light_mat/matexpr/mat_arith.h>
#include <light_mat/matexpr/mat_emath.h>
#include <light_mat/matexpr/mat_pred.h>
#include <light_mat/matexpr/mat_cast.h>
#include <light_mat/matexpr/repvec_expr.h>

#include <light_mat/mateval/mat_reduce.h>
#include <light_mat/mateval/mat_enorms.h>
#include <light_mat/mateval/mat_allany.h>

namespace dolphin
{
	using lmat::type_;

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
	using lmat::ref_row;
	using lmat::ref_grid;

	using lmat::cref_matrix;
	using lmat::cref_col;
	using lmat::cref_row;
	using lmat::cref_grid;

	using lmat::ref_block;
	using lmat::ref_grid;
	using lmat::step_col;
	using lmat::step_row;

	using lmat::cref_block;
	using lmat::cref_grid;
	using lmat::cstep_col;
	using lmat::cstep_row;

	using lmat::meta::is_contiguous;
	using lmat::meta::is_percol_contiguous;
	using lmat::meta::supports_linear_index;

	using lmat::is_vector;
	using lmat::is_column;
	using lmat::is_row;
	using lmat::is_empty;
	using lmat::is_scalar;
	using lmat::is_square;

	using lmat::as_col;
	using lmat::as_row;

	using lmat::in_;
	using lmat::out_;
	using lmat::in_out_;

	using lmat::ewise;
	using lmat::map;
	using lmat::map_to;
	using lmat::accum_to;

	// import math functions

	namespace math = lmat::math;

	using lmat::repcol;
	using lmat::reprow;

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

	// import reduction functions

	using lmat::fold;

	using lmat::sum;
	using lmat::mean;
	using lmat::maximum;
	using lmat::minimum;

	using lmat::sqsum;
	using lmat::asum;
	using lmat::amax;

	using lmat::colwise_sum;
	using lmat::colwise_mean;
	using lmat::colwise_maximum;
	using lmat::colwise_minimum;

	using lmat::colwise_sqsum;
	using lmat::colwise_asum;
	using lmat::colwise_amax;

	using lmat::rowwise_sum;
	using lmat::rowwise_mean;
	using lmat::rowwise_maximum;
	using lmat::rowwise_minimum;

	using lmat::rowwise_sqsum;
	using lmat::rowwise_asum;
	using lmat::rowwise_amax;

	namespace norms = lmat::norms;
	using lmat::norm;
	using lmat::colwise_norm;
	using lmat::rowwise_norm;

	using lmat::all;
	using lmat::any;
	using lmat::colwise_all;
	using lmat::colwise_any;

}

#endif
