//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#pragma once

#include <El.hpp>

#include <cassert>
#include <iostream>
#include <vector>

// FIXME: Use boost::math::tools::polynomial instead

// A univariate polynomial
//
//   p(x) = a_0 + a_1 x + a_2 x^2 + ... + a_n x^n
//
class Polynomial
{
public:
  // Coefficients {a_0, a_1, ..., a_n} in increasing order of degree
  std::vector<El::BigFloat> coefficients;

  // The zero polynomial
  Polynomial() : coefficients(1, 0) {}
  Polynomial(const size_t &size, const El::BigFloat &default_element)
      : coefficients(size, default_element)
  {}

  // Degree of p(x)
  int64_t degree() const { return coefficients.size() - 1; };

  // Evaluate p(x) for some x using horner's method
  El::BigFloat operator()(const El::BigFloat &x) const
  {
    assert(!coefficients.empty());
    auto coefficient(coefficients.rbegin());
    El::BigFloat result(*coefficient);
    ++coefficient;
    for(; coefficient != coefficients.rend(); ++coefficient)
      {
        result *= x;
        result += *coefficient;
      }
    return result;
  }

  // Print p(x), for debugging purposes
  friend std::ostream &operator<<(std::ostream &os, const Polynomial &p)
  {
    for(int i = p.degree(); i >= 0; i--)
      {
        os << p.coefficients[i];
        if(i > 1)
          {
            os << "x^" << i << " + ";
          }
        else if(i == 1)
          {
            os << "x + ";
          }
      }
    return os;
  }

  inline void operator+=(const El::BigFloat &b)
  {
	  coefficients[0] += b;
  }
  inline void operator-=(const El::BigFloat &b)
  {
	  coefficients[0] -= b;
  }
  inline void operator*=(const El::BigFloat &b)
  {
	  for (auto &coeff : coefficients) coeff *= b;
  }
  inline void operator/=(const El::BigFloat &b)
  {
	  for (auto &coeff : coefficients) coeff /= b;
  }
  inline void operator+=(const Polynomial &b)
  {
	  if(coefficients.size() < b.coefficients.size())
		  coefficients.resize(b.coefficients.size(), 0);
	  for (int i = 0; i < b.coefficients.size(); i++)
		  coefficients[i] += b.coefficients[i];
  }
  inline void operator-=(const Polynomial &b)
  {
	  if (coefficients.size() < b.coefficients.size())
		  coefficients.resize(b.coefficients.size(), 0);
	  for (int i = 0; i < b.coefficients.size(); i++)
		  coefficients[i] -= b.coefficients[i];
  }
  inline void operator-()
  {
	  for (auto &coefficient : coefficients)
		  coefficient = -coefficient;
  }

  // follow https://www.dealii.org/current/doxygen/deal.II/polynomial_8cc_source.html#l00572
  void shift(const El::BigFloat offset)
  {
	  // Copy coefficients to a vector of
	  // accuracy given by the argument
	  std::vector<El::BigFloat> new_coefficients(coefficients.begin(), coefficients.end());

	  El::BigFloat f_RHS;

	  // Traverse all coefficients from
	  // c_1. c_0 will be modified by
	  // higher degrees, only.
	  for (unsigned int d = 1; d < new_coefficients.size(); ++d)
	  {
		  const unsigned int n = d;
		  // Binomial coefficients are
		  // needed for the
		  // computation. The rightmost
		  // value is unity.
		  mpz_class binomial_coefficient = 1;

		  // Powers of the offset will be
		  // needed and computed
		  // successively.
		  El::BigFloat offset_power = offset;

		  // Compute (x+offset)^d
		  // and modify all values c_k
		  // with k<d.
		  // The coefficient in front of
		  // x^d is not modified in this step.
		  for (unsigned int k = 0; k < d; ++k)
		  {
			  // Recursion from Bronstein
			  // Make sure no remainders
			  // occur in integer
			  // division.
			  binomial_coefficient = (binomial_coefficient * (n - k)) / (k + 1);

			  f_RHS.gmp_float = new_coefficients[d].gmp_float * offset_power.gmp_float * binomial_coefficient;
			  // f_RHS = new_coefficients[d] * binomial_coefficient * offset_power;

			  new_coefficients[d - k - 1] += f_RHS;
			  offset_power *= offset;
		  }
		  // The binomial coefficient
		  // should have gone through a
		  // whole row of Pascal's
		  // triangle.
		  if (binomial_coefficient != 1) 
		  {
			  std::stringstream ss;
			  ss << "Polynomial::shift internal error. binomial_coefficient=" << binomial_coefficient 
				  << ", n=" << n << ", d=" << d << "\n";
			  throw std::runtime_error(ss.str());
		  }
	  }

	  // copy new elements to old vector
	  coefficients.assign(new_coefficients.begin(), new_coefficients.end());
  }


};

inline Polynomial operator/(const Polynomial &a, const El::BigFloat &b)
{
  Polynomial result(0, 0);
  result.coefficients.reserve(a.degree() + 1);
  const El::BigFloat inverse(1 / b);
  for(auto &coefficient : a.coefficients)
    {
      result.coefficients.emplace_back(coefficient * inverse);
    }
  return result;
}


// Convenience functions to avoid copies
inline void swap(std::vector<Polynomial> &polynomials,
                 std::vector<std::vector<El::BigFloat>> &elements_vector)
{
  for(auto &elements : elements_vector)
    {
      polynomials.emplace_back();
      std::swap(polynomials.back().coefficients, elements);
    }
}

inline void swap(
  std::vector<std::vector<Polynomial>> &polynomials_vector,
  std::vector<std::vector<std::vector<El::BigFloat>>> &elements_vector_vector)
{
  for(auto &elements_vector : elements_vector_vector)
    {
      polynomials_vector.emplace_back();
      swap(polynomials_vector.back(), elements_vector);
    }
}

inline void swap(
  std::vector<std::vector<std::vector<Polynomial>>> &polynomials_vector_vector,
  std::vector<std::vector<std::vector<std::vector<El::BigFloat>>>>
    &elements_vector_vector_vector)
{
  for(auto &elements_vector_vector : elements_vector_vector_vector)
    {
      polynomials_vector_vector.emplace_back();
      swap(polynomials_vector_vector.back(), elements_vector_vector);
    }
}
