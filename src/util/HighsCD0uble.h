/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file util/HighsCD0uble.h
 * @brief Quad precision type implemented with two standard double precision
 *        representing the value and a compensation term
 */
#ifndef UTIL_HIGHSCD0UBLE_H_
#define UTIL_HIGHSCD0UBLE_H_

#include <cmath>
#include <cstdint>

/// A compensated double number achieving roughly quad precision on the
/// supported operations

class HighsCD0uble {
 private:
  double hi;
  double lo;

  // The following functions are implemented as described in:
  // Rump, Siegfried M. "High precision evaluation of nonlinear functions."
  // Proceedings of. 2005.

  /// performs an exact transformation such that x + y = a + b
  /// and x = double(a + b). The operation uses 6 flops (addition/substraction).
  static void two_sum(double& x, double& y, double a, double b) {
    x = a + b;
    double z = x - a;
    y = (a - (x - z)) + (b - z);
  }

  /// splits a 53 bit double precision number into two 26 bit parts
  /// such that x + y = a holds exactly
  static void split(double& x, double& y, double a) {
    constexpr double factor = double((1 << 27) + 1);
    double c = factor * a;
    x = c - (c - a);
    y = a - x;
  }

  /// performs an exact transformation such that x + y = a * b
  /// and x = double(a * b). The operation uses 10 flops for
  /// addition/substraction and 7 flops for multiplication.
  static void two_product(double& x, double& y, double a, double b) {
    x = a * b;
    double a1, a2, b1, b2;
    split(a1, a2, a);
    split(b1, b2, b);
    y = a2 * b2 - (((x - a1 * b1) - a2 * b1) - a1 * b2);
  }

  HighsCD0uble(double hi, double lo) : hi(hi), lo(lo) {}

 public:
  HighsCD0uble() = default;

  HighsCD0uble(double val) : hi(val), lo(0.0) {}

  explicit operator double() const { return hi + lo; }

  HighsCD0uble& operator+=(double v) {
    double c;
    two_sum(hi, c, v, hi);
    lo += c;
    return *this;
  }

  HighsCD0uble& operator+=(const HighsCD0uble& v) {
    (*this) += v.hi;
    lo += v.lo;
    return *this;
  }

  HighsCD0uble& operator-=(double v) {
    (*this) += -v;
    return *this;
  }

  HighsCD0uble& operator-=(const HighsCD0uble& v) {
    (*this) -= v.hi;
    lo -= v.lo;
    return *this;
  }

  HighsCD0uble& operator*=(double v) {
    double c = lo * v;
    two_product(hi, lo, hi, v);
    *this += c;
    return *this;
  }

  HighsCD0uble& operator*=(const HighsCD0uble& v) {
    double c1 = hi * v.lo;
    double c2 = lo * v.hi;
    two_product(hi, lo, hi, v.hi);
    *this += c1;
    *this += c2;
    return *this;
  }

  HighsCD0uble& operator/=(double v) {
    HighsCD0uble d(hi / v, lo / v);
    HighsCD0uble c = d * v - (*this);
    c.hi /= v;
    c.lo /= v;
    *this = d - c;
    return *this;
  }

  HighsCD0uble& operator/=(const HighsCD0uble& v) {
    double vdbl = v.hi + v.lo;
    HighsCD0uble d(hi / vdbl, lo / vdbl);
    HighsCD0uble c = d * v - (*this);
    c.hi /= vdbl;
    c.lo /= vdbl;
    *this = d - c;
    return *this;
  }

  HighsCD0uble operator-() const { return HighsCD0uble(-hi, -lo); }

  HighsCD0uble operator+(double v) const {
    HighsCD0uble res;

    two_sum(res.hi, res.lo, hi, v);
    res.lo += lo;

    return res;
  }

  HighsCD0uble operator+(const HighsCD0uble& v) const {
    HighsCD0uble res = (*this) + v.hi;
    res.lo += v.lo;

    return res;
  }

  friend HighsCD0uble operator+(double a, const HighsCD0uble& b) {
    return b + a;
  }

  HighsCD0uble operator-(double v) const {
    HighsCD0uble res;

    two_sum(res.hi, res.lo, hi, -v);
    res.lo += lo;

    return res;
  }

  HighsCD0uble operator-(const HighsCD0uble& v) const {
    HighsCD0uble res = (*this) - v.hi;
    res.lo -= v.lo;

    return res;
  }

  friend HighsCD0uble operator-(double a, const HighsCD0uble& b) {
    return -b + a;
  }

  HighsCD0uble operator*(double v) const {
    HighsCD0uble res;

    two_product(res.hi, res.lo, hi, v);
    res += lo * v;

    return res;
  }

  HighsCD0uble operator*(const HighsCD0uble& v) const {
    HighsCD0uble res = (*this) * v.hi;
    res += hi * v.lo;

    return res;
  }

  friend HighsCD0uble operator*(double a, const HighsCD0uble& b) {
    return b * a;
  }

  HighsCD0uble operator/(double v) const {
    HighsCD0uble res = *this;

    res /= v;

    return res;
  }

  HighsCD0uble operator/(const HighsCD0uble& v) const {
    HighsCD0uble res = (*this);

    res /= v;

    return res;
  }

  friend HighsCD0uble operator/(double a, const HighsCD0uble& b) {
    return HighsCD0uble(a) / b;
  }

  bool operator>(const HighsCD0uble& other) const {
    return double(*this) > double(other);
  }

  bool operator>(double other) const { return double(*this) > other; }

  friend bool operator>(double a, const HighsCD0uble& b) {
    return a > double(b);
  }

  bool operator<(const HighsCD0uble& other) const {
    return double(*this) < double(other);
  }

  bool operator<(double other) const { return double(*this) < other; }

  friend bool operator<(double a, const HighsCD0uble& b) {
    return a < double(b);
  }

  bool operator>=(const HighsCD0uble& other) const {
    return double(*this) >= double(other);
  }

  bool operator>=(double other) const { return double(*this) >= other; }

  friend bool operator>=(double a, const HighsCD0uble& b) {
    return a >= double(b);
  }

  bool operator<=(const HighsCD0uble& other) const {
    return double(*this) <= double(other);
  }

  bool operator<=(double other) const { return double(*this) <= other; }

  friend bool operator<=(double a, const HighsCD0uble& b) {
    return a <= double(b);
  }

  bool operator==(const HighsCD0uble& other) const {
    return double(*this) == double(other);
  }

  bool operator==(double other) const { return double(*this) == other; }

  friend bool operator==(double a, const HighsCD0uble& b) {
    return a == double(b);
  }

  bool operator!=(const HighsCD0uble& other) const {
    return double(*this) != double(other);
  }

  bool operator!=(double other) const { return double(*this) != other; }

  friend bool operator!=(double a, const HighsCD0uble& b) {
    return a != double(b);
  }

  void renormalize() { two_sum(hi, lo, hi, lo); }

  friend HighsCD0uble abs(const HighsCD0uble& v) { return v < 0 ? -v : v; }
  friend HighsCD0uble fabs(const HighsCD0uble& v) { return v < 0 ? -v : v; }

  friend HighsCD0uble sqrt(const HighsCD0uble& v) {
    double c = std::sqrt(v.hi + v.lo);

    // guard against division by zero
    if (c == 0.0) return 0.0;

    // calculate precise square root by newton step
    HighsCD0uble res = v / c;
    res += c;
    // multiplication by 0.5 is exact
    res.hi *= 0.5;
    res.lo *= 0.5;
    return res;
  }

  friend HighsCD0uble floor(const HighsCD0uble& x) {
    double floor_x = std::floor(double(x));
    HighsCD0uble res;

    two_sum(res.hi, res.lo, floor_x, std::floor(double(x - floor_x)));
    return res;
  }

  friend HighsCD0uble ceil(const HighsCD0uble& x) {
    double ceil_x = std::ceil(double(x));
    HighsCD0uble res;

    two_sum(res.hi, res.lo, ceil_x, std::ceil(double(x - ceil_x)));
    return res;
  }

  friend HighsCD0uble round(const HighsCD0uble& x) { return floor(x + 0.5); }

};

#endif
