/**
 * @file properties.h
 *
 * @brief Emulation of property syntax
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef DOLPHIN_PROPERTIES_H_
#define DOLPHIN_PROPERTIES_H_

#include <dolphin/common/common_base.h>
#include <functional>

namespace dolphin
{

	class invalid_property_value : public invalid_argument
	{
	public:
		invalid_property_value(const char *msg)
		: invalid_argument(msg) { }
	};


	template<typename T, class Validator>
	DOLPHIN_ENSURE_INLINE
	inline const T& validated_property_value(const T& v, const Validator& validator, const char *errmsg)
	{
		if (!validator(v))
			throw invalid_property_value(errmsg);

		return v;
	}


	template<typename T>
	struct trivial_validator
	{
		DOLPHIN_ENSURE_INLINE
		bool operator() (const T& v) const
		{
			return true;
		}
	};


	template<typename T>
	class simple_property
	{
	public:
		typedef T value_type;

		explicit simple_property(const T& init)
		: m_value(init), m_validator(trivial_validator<T>()), m_errmsg("")
		{ }

		template<class Validator>
		explicit simple_property(const T& init,
				const Validator& validator, const char *errmsg)
		: m_value(validated_property_value(init, validator, errmsg))
		, m_validator(validator), m_errmsg(errmsg) { }

		DOLPHIN_ENSURE_INLINE
		operator T() const
		{
			return m_value;
		}

		DOLPHIN_ENSURE_INLINE
		const T& get() const
		{
			return m_value;
		}

		bool validate(const T& v) const
		{
			return m_validator(v);
		}

		void set(const T& v)
		{
			if (!validate(v))
				throw invalid_property_value(m_errmsg);
			m_value = v;
		}

	private:
		T m_value;
		std::function<bool(const T&)> m_validator;
		const char *m_errmsg;
	};


	/********************************************
	 *
	 *  commonly useful validators
	 *
	 ********************************************/

	template<typename T>
	struct require_gt_
	{
		const T bound;

		explicit require_gt_(const T& b)
		: bound(b) { }

		DOLPHIN_ENSURE_INLINE
		bool operator() (const T& v) const
		{
			return v > bound;
		}
	};

	template<typename T>
	struct require_ge_
	{
		const T bound;

		DOLPHIN_ENSURE_INLINE
		explicit require_ge_(const T& b)
		: bound(b) { }

		DOLPHIN_ENSURE_INLINE
		bool operator() (const T& v) const
		{
			return v >= bound;
		}
	};

	template<typename T>
	struct require_lt_
	{
		const T bound;

		DOLPHIN_ENSURE_INLINE
		explicit require_lt_(const T& b)
		: bound(b) { }

		DOLPHIN_ENSURE_INLINE
		bool operator() (const T& v) const
		{
			return v < bound;
		}
	};

	template<typename T>
	struct require_le_
	{
		const T bound;

		DOLPHIN_ENSURE_INLINE
		explicit require_le_(const T& b)
		: bound(b) { }

		DOLPHIN_ENSURE_INLINE
		bool operator() (const T& v) const
		{
			return v <= bound;
		}
	};

	template<typename T>
	struct require_within_
	{
		const T lbound;
		const T ubound;

		DOLPHIN_ENSURE_INLINE
		require_within_(const T& lb, const T& ub)
		: lbound(lb), ubound(ub) { }

		DOLPHIN_ENSURE_INLINE
		bool operator() (const T& v) const
		{
			return v >= lbound && v <= ubound;
		}
	};


	template<typename T>
	DOLPHIN_ENSURE_INLINE
	require_gt_<T> require_gt(const T& b)
	{
		return require_gt_<T>(b);
	}

	template<typename T>
	DOLPHIN_ENSURE_INLINE
	require_ge_<T> require_ge(const T& b)
	{
		return require_ge_<T>(b);
	}

	template<typename T>
	DOLPHIN_ENSURE_INLINE
	require_lt_<T> require_lt(const T& b)
	{
		return require_lt_<T>(b);
	}

	template<typename T>
	DOLPHIN_ENSURE_INLINE
	require_le_<T> require_le(const T& b)
	{
		return require_le_<T>(b);
	}

	template<typename T>
	DOLPHIN_ENSURE_INLINE
	require_within_<T> require_within(const T& a, const T& b)
	{
		return require_within_<T>(a, b);
	}

}

#endif
