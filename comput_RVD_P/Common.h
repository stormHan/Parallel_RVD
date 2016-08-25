/*
common definition
*/
#ifndef H_COMMON_H
#define H_COMMON_H

namespace P_RVD
{
	/*
	the type for manipulating indeces
	*/
	typedef unsigned int t_index;

	typedef int signed_t_index;

	/*
	for glut
	*/
	typedef unsigned int uint;

	/*
	find the min value
	*/
	template <class T>
	T inline geo_min(T& a, T& b)
	{
		return (a < b) ? a : b;
	}

	/*
	find the max value
	*/
	template <class T>
	T inline geo_max(T& a, T& b)
	{
		return (a < b) ? b : a;
	}

	/*
	swap function
	*/
	template <class T>
	inline void geo_swap(T& a, T& b)
	{
		T c = a;
		a = b; b = c;
	}

	/**
	* \brief Integer constants that represent the sign of a value
	*/
	enum Sign {
		/** Value is negative */
		NEGATIVE = -1,
		/** Value is zero */
		ZERO = 0,
		/** Value is positive */
		POSITIVE = 1
	};

	/**
	* \brief Gets the sign of a value
	* \details Returns -1, 0, or 1 whether value \p x is resp. negative, zero
	* or positive. The function uses operator<() and operator>() to compare
	* the value agains to 0 (zero). The integer constant zero must make
	* senses for the type of the value, or T must be constructible from
	* integer constant zero.
	* \param[in] x the value to test
	* \tparam T the type of the value
	* \return the sign of the value
	* \see Sign
	*/
	template <class T>
	inline Sign geo_sgn(const T& x) {
		return (x > 0) ? POSITIVE : (
			(x < 0) ? NEGATIVE : ZERO
			);
	}
}

#endif /* H_COMMON_H */