#pragma once

#include <cstddef>
#include <stdexcept>
#include <sstream>

// Extracted from Boost
// boost/math/special_functions/chebyshev.hpp

//  (C) Copyright Nick Thompson 2017.
//  Use, modification and distribution are subject to the Boost
//  Software License, Version 1.0.

// Boost Software License - Version 1.0 - August 17th, 2003

// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:

// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

namespace detail {
template<class Real>
inline Real unchecked_chebyshev_clenshaw_recurrence(const Real* const c, size_t length, const Real & a, const Real & b, const Real& x)
{
    Real t;
    Real u;
    // This cutoff is not super well defined, but it's a good estimate.
    // See "An Error Analysis of the Modified Clenshaw Method for Evaluating Chebyshev and Fourier Series"
    // J. OLIVER, IMA Journal of Applied Mathematics, Volume 20, Issue 3, November 1977, Pages 379–391
    // https://doi.org/10.1093/imamat/20.3.379
    const Real cutoff = 0.6;
    if (x - a < b - x)
    {
        u = 2*(x-a)/(b-a);
        t = u - 1;
        if (t > -cutoff)
        {
            Real b2 = 0;
            Real b1 = c[length -1];
            for(size_t j = length - 2; j >= 1; --j)
            {
                Real tmp = 2*t*b1 - b2 + c[j];
                b2 = b1;
                b1 = tmp;
            }
            return t*b1 - b2 + c[0]/2;
        }
        else
        {
            Real b = c[length -1];
            Real d = b;
            Real b2 = 0;
            for (size_t r = length - 2; r >= 1; --r)
            {
                d = 2*u*b - d + c[r];
                b2 = b;
                b = d - b;
            }
            return t*b - b2 + c[0]/2;
        }
    }
    else
    {
        u = -2*(b-x)/(b-a);
        t = u + 1;
        if (t < cutoff)
        {
            Real b2 = 0;
            Real b1 = c[length -1];
            for(size_t j = length - 2; j >= 1; --j)
            {
                Real tmp = 2*t*b1 - b2 + c[j];
                b2 = b1;
                b1 = tmp;
            }
            return t*b1 - b2 + c[0]/2;
        }
        else
        {
            Real b = c[length -1];
            Real d = b;
            Real b2 = 0;
            for (size_t r = length - 2; r >= 1; --r)
            {
                d = 2*u*b + d + c[r];
                b2 = b;
                b = d + b;
            }
            return t*b - b2 + c[0]/2;
        }
    }
}

} // namespace detail

template<class Real>
inline Real chebyshev_clenshaw_recurrence(const Real* const c, size_t length, const Real & a, const Real & b, const Real& x)
{
    if (x < a || x > b)
    {
      std::stringstream ss;
      ss << "x in [a, b] is required: x = "
         << x
         << ", a = "
         << a
         << ", b = "
         << b;
      throw std::domain_error(ss.str());
    }
    if (length < 2)
    {
        if (length == 0)
        {
            return 0;
        }
        return c[0]/2;
    }
    return detail::unchecked_chebyshev_clenshaw_recurrence(c, length, a, b, x);
}