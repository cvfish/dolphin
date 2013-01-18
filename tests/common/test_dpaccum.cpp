/**
 * @file test_dpaccum.cpp
 *
 * Unit testing of dispatched accumulation
 * 
 * @author Dahua Lin 
 */

#include "../test_base.h"
#include <dolphin/common/dpaccum.h>

using namespace dolphin;
using namespace dolphin::test;

SIMPLE_CASE( test_add_counts )
{
	const index_t K = 12;
	const index_t len = 200;
	dense_col<uint32_t> cnts(K, zero());
	dense_col<index_t> I(len, zero());

	fill_randi(I, (index_t)0, K+2);
	add_counts(I, cnts);

	dense_col<uint32_t> c0(K, zero());
	for (index_t i = 0; i < len; ++i)
	{
		index_t k = I[i];
		if (k >= 0 && k < K) c0[k]++;
	}

	ASSERT_VEC_EQ(K, cnts, c0);
}

AUTO_TPACK( test_counts )
{
	ADD_SIMPLE_CASE( test_add_counts )
}

