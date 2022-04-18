#define QDOUBLE_IMPLEMENTATION
#include "qdouble.h"

#include "test_assert.h"

#include <algorithm>
#include <iostream>

namespace
  {

  inline double get_epsilon()
    {
    double eps = 1.0;
    while (1.0 != (eps + 1.0))
      eps /= 2.0;
    return eps * 2.0;
    }

  void test_two_sum_1()
    {
    double err;
    double s = two_sum(1.0, 2.0, err);
    TEST_EQ(3.0, s);
    TEST_EQ(0.0, err);
    }

  void test_two_sum_2()
    {
    double err;
    double s = two_sum(1.0, 1e-16, err);
    TEST_EQ(1.0, s);
    TEST_EQ(1e-16, err);
    }

  void test_two_sum_3()
    {
    double eps_div_2 = get_epsilon() / 2.0;
    double err;
    double s = two_sum(1.0, eps_div_2, err);
    TEST_EQ(1.0, s);
    TEST_EQ(eps_div_2, err);
    }

  void test_two_sum_4()
    {
    double eps = get_epsilon();
    double err;
    double s = two_sum(1.0, eps, err);
    TEST_EQ(1.0 + eps, s);
    TEST_EQ(0.0, err);
    }
  
  void test_simple_arithmetic()
    {
    qdouble a(1.0);
    qdouble b(2.0);

    double c(3.0);

    TEST_EQ(4.0, a + c);
    TEST_EQ(4.0, c + a);

    TEST_EQ(3.0, a + b);
    TEST_EQ(3.0, b + a);


    TEST_EQ(-2.0, a - c);
    TEST_EQ(2.0, c - a);

    TEST_EQ(-1.0, a - b);
    TEST_EQ(1.0, b - a);
    }

  void test_output_test()
    {
    std::cout.precision(60);
    std::cout << qdouble_pi << std::endl;
    qdouble x = make_qdouble("1.0");
    x /= 3.0;
    qdouble y;
    y = pow(qdouble(2.0), 3);
    std::cout << "y = " << y << std::endl;
    std::cout << "x = " << x << std::endl;


    qdouble a;
    qdouble b = make_qdouble("0.1");
    std::cout << b << std::endl;
    a = sqrt(b);
    std::cout << " sqrt(0.1) = " << a << std::endl;
    std::cout << " sqrt(0.1) * sqrt(0.1) = " << a * a << std::endl;
    }

  void test_polynomial()
    {
    std::cout << std::endl;
    std::cout << "Test 1.  (Polynomial)." << std::endl;

    int n = 8;
    qdouble* c = new qdouble[n];
    qdouble x, y;

    for (int i = 0; i < n; ++i)
      c[i] = static_cast<double>(i + 1);

    x = polyroot(c, n - 1, qdouble(0.0));
    y = polyeval(c, n - 1, x);

    std::cout.precision(6);
    std::cout << "Root Found:  x  = " << x << std::endl;
    std::cout << "           p(x) = " << y << std::endl;

    delete[] c;
    TEST_ASSERT(to_double(y) < 4.0 * qdouble_eps);
    }

  void test_machin_formula_pi()
    {
    std::cout << std::endl;
    std::cout << "Test 2.  (Machin's Formula for Pi)." << std::endl;

    // Use the Machin's arctangent formula:
    //
    //pi / 4  =  4 arctan(1/5) - arctan(1/239)
    //
    //The arctangent is computed based on the Taylor series expansion
    //
    //arctan(x) = x - x^3 / 3 + x^5 / 5 - x^7 / 7 + ...
    //

    qdouble s1, s2, t, r;
    int k;
    int sign;
    double d;
    double err;

    // Compute arctan(1/5)
    d = 1.0;
    t = qdouble(1.0) / 5.0;
    r = sqr(t);
    s1 = 0.0;
    k = 0;

    sign = 1;
    while (t > qdouble_eps)
      {
      ++k;
      if (sign < 0)
        s1 -= (t / d);
      else
        s1 += (t / d);

      d += 2.0;
      t *= r;
      sign = -sign;
      }

    std::cout << k << " Iterations" << std::endl;

    // Compute arctan(1/239) 
    d = 1.0;
    t = qdouble(1.0) / 239.0;
    r = sqr(t);
    s2 = 0.0;
    k = 0;

    sign = 1;
    while (t > qdouble_eps)
      {
      ++k;
      if (sign < 0)
        s2 -= (t / d);
      else
        s2 += (t / d);

      d += 2.0;
      t *= r;
      sign = -sign;
      }

    std::cout << k << " Iterations" << std::endl;

    qdouble p = 4.0 * s1 - s2;

    p *= 4.0;
    err = abs(to_double(p - qdouble_pi));


    std::cout.precision(qdouble_ndigits);
    std::cout << "   pi = " << p << std::endl;
    std::cout << "  _pi = " << qdouble_pi << std::endl;

    std::cout.precision(6);
    std::cout << "error = " << err << " = " << err / qdouble_eps << " eps" << std::endl;


    TEST_ASSERT(err < 8.0 * qdouble_eps);
    }

  void test_salamin_brent_quadratic_formula_pi()
    {
    std::cout << std::endl;
    std::cout << "Test 3.  (Salamin-Brent Quadratic Formula for Pi)." << std::endl;
    std::cout.precision(qdouble_ndigits);

    qdouble a, b, s, p;
    qdouble a_new, b_new, p_old;
    double m;
    double err;
    const int max_iter = 20;

    a = 1.0;
    b = sqrt(qdouble(0.5));
    s = 0.5;
    m = 1.0;

    p = 2.0 * sqr(a) / s;
    std::cout << "Iteration  0: " << p << std::endl;
    for (int i = 1; i <= max_iter; ++i)
      {
      m *= 2.0;
      a_new = 0.5 * (a + b);
      b_new = a * b;
      s -= m * (sqr(a_new) - b_new);
      a = a_new;
      b = sqrt(b_new);
      p_old = p;
      p = 2.0 * sqr(a) / s;
      std::cout << "Iteration " << std::setw(2) << i << ": " << p << std::endl;
      if (abs(to_double(p - p_old)) < 64 * qdouble_eps)
        break;
      }

    err = abs(to_double(p - qdouble_pi));

    std::cout << "         _pi: " << qdouble_pi << std::endl;
    std::cout.precision(6);
    std::cout << "       error: " << err << " = " << err / qdouble_eps << " eps" << std::endl;


    // for some reason, this test gives relatively large error compared
    // to other tests.  May need to be looked at more closely.
    TEST_ASSERT(err < 1024.0 * qdouble_eps);
    }

  void test_borwein_quartic_formula_for_pi()
    {
    std::cout << std::endl;
    std::cout << "Test 4.  (Borwein Quartic Formula for Pi)." << std::endl;
    std::cout.precision(qdouble_ndigits);

    qdouble a, y, p, r, p_old;
    double m;
    double err;
    const int max_iter = 20;

    a = 6.0 - 4.0 * sqrt(qdouble(2.0));
    y = sqrt(qdouble(2.0)) - 1.0;
    m = 2.0;

    p = 1.0 / a;

    std::cout << "Iteration  0: " << p << std::endl;

    for (int i = 1; i <= max_iter; ++i)
      {
      m *= 4.0;
      r = nroot(1.0 - sqr(sqr(y)), 4);
      y = (1.0 - r) / (1.0 + r);
      a = a * sqr(sqr(1.0 + y)) - m * y * (1.0 + y + sqr(y));

      p_old = p;
      p = 1.0 / a;
      std::cout << "Iteration " << std::setw(2) << i << ": " << p << std::endl;
      if (abs(to_double(p - p_old)) < 16 * qdouble_eps)
        break;
      }

    err = abs(to_double(p - qdouble_pi));

    std::cout << "         _pi: " << qdouble_pi << std::endl;
    std::cout.precision(6);
    std::cout << "       error: " << err << " = " << err / qdouble_eps << " eps" << std::endl;


    TEST_ASSERT(err < 256.0 * qdouble_eps);
    }

  void test_taylor_series_for_e()
    {
    std::cout << std::endl;
    std::cout << "Test 5.  (Taylor Series Formula for E)." << std::endl;
    std::cout.precision(qdouble_ndigits);

    // Use Taylor series
    //
    // e = 1 + 1 + 1/2! + 1/3! + 1/4! + ...
    //
    // To compute e.
    //

    qdouble s = 2.0, t = 1.0;
    double n = 1.0;
    double delta;
    int i = 0;

    while (t > qdouble_eps)
      {
      ++i;
      n += 1.0;
      t /= n;
      s += t;
      }

    delta = abs(to_double(s - qdouble_e));


    std::cout << "    e = " << s << std::endl;
    std::cout << "   _e = " << qdouble_e << std::endl;

    std::cout.precision(6);
    std::cout << "error = " << delta << " = " << delta / qdouble_eps << " eps" << std::endl;
    std::cout << i << " iterations." << std::endl;

    TEST_ASSERT(delta < 64.0 * qdouble_eps);
    }

  void test_taylor_series_for_log_2()
    {
    std::cout << std::endl;
    std::cout << "Test 6.  (Taylor Series Formula for Log 2)." << std::endl;
    std::cout.precision(qdouble_ndigits);

    // Use the Taylor series
    //
    //-log(1-x) = x + x^2/2 + x^3/3 + x^4/4 + ...
    //
    //with x = 1/2 to get  log(1/2) = -log 2.
    //

    qdouble s = 0.5;
    qdouble t = 0.5;
    double delta;
    double n = 1.0;
    double i = 0;

    while (abs(t) > qdouble_eps)
      {
      ++i;
      n += 1.0;
      t *= 0.5;
      s += (t / n);
      }

    delta = abs(to_double(s - qdouble_log2));

    std::cout << " log2 = " << s << std::endl;
    std::cout << "_log2 = " << qdouble_log2 << std::endl;
    std::cout.precision(6);
    std::cout << "error = " << delta << " = " << (delta / qdouble_eps) << " eps" << std::endl;
    std::cout << i << " iterations." << std::endl;


    TEST_ASSERT(delta < 4.0 * qdouble_eps);
    }

  void test_sanity_check_exp()
    {
    std::cout << std::endl;
    std::cout << "Test 7.  (Sanity TEST_ASSERT for exp)." << std::endl;
    std::cout.precision(qdouble_ndigits);

    // Do simple sanity TEST_ASSERT
    //
    //  e^2 = exp(2)
    //      = exp(-13/4) * exp(-9/4) * exp(-5/4) * exp(-1/4) *
    //        exp(3/4) * exp(7/4) * exp(11/4) * exp(15/4)
    //

    qdouble t = -3.25;
    qdouble p = 1.0;

    for (int i = 0; i < 8; i++, t += 1.0)
      {
      p = p * exp(t);
      }

    qdouble t1 = exp(qdouble(2.0));
    qdouble t2 = sqr(qdouble_e);
    double delta = std::max(abs(to_double(t1 - p)), abs(to_double(t2 - p)));


    std::cout << "result = " << p << std::endl;
    std::cout << "exp(2) = " << t1 << std::endl;
    std::cout << "   e^2 = " << t2 << std::endl;
    std::cout.precision(6);
    std::cout << " error = " << delta << " = " << (delta / qdouble_eps) << " eps" << std::endl;

    TEST_ASSERT(delta < 16.0 * qdouble_eps);
    }

  void test_sanity_check_sin_cos()
    {
    std::cout << std::endl;
    std::cout << "Test 8.  (Sanity TEST_ASSERT for sin / cos)." << std::endl;
    std::cout.precision(qdouble_ndigits);

    // Do simple sanity TEST_ASSERT
    //
    // sin(x) = sin(5x/7)cos(2x/7) + cos(5x/7)sin(2x/7)
    //
    // cos(x) = cos(5x/7)cos(2x/7) - sin(5x/7)sin(2x/7);
    //

    qdouble x = qdouble_pi / 3.0;
    qdouble x1 = 5.0 * x / 7.0;
    qdouble x2 = 2.0 * x / 7.0;

    qdouble r1 = sin(x1) * cos(x2) + cos(x1) * sin(x2);
    qdouble r2 = cos(x1) * cos(x2) - sin(x1) * sin(x2);
    qdouble t1 = sqrt(qdouble(3.0)) / 2.0;
    qdouble t2 = 0.5;

    double delta = std::max(abs(to_double(t1 - r1)), abs(to_double(t2 - r2)));


    std::cout << "  r1 = " << r1 << std::endl;
    std::cout << "  t1 = " << t1 << std::endl;
    std::cout << "  r2 = " << r2 << std::endl;
    std::cout << "  t2 = " << t2 << std::endl;
    std::cout.precision(6);
    std::cout << " error = " << delta << " = " << (delta / qdouble_eps) << " eps" << std::endl;


    TEST_ASSERT(delta < 4.0 * qdouble_eps);
    }

  void test_make_qdouble()
    {
    qdouble q = make_qdouble("12e15");
    TEST_EQ(12e15, q.a[0]);
    TEST_EQ(0.0, q.a[1]);
    TEST_EQ(0.0, q.a[2]);
    TEST_EQ(0.0, q.a[3]);
    q = make_qdouble("12edf");
    TEST_ASSERT(q != q); // not a number
    q = make_qdouble("12e0");
    TEST_EQ(12, q.a[0]);
    TEST_EQ(0.0, q.a[1]);
    TEST_EQ(0.0, q.a[2]);
    TEST_EQ(0.0, q.a[3]);
    q = make_qdouble("12e0000");
    TEST_EQ(12, q.a[0]);
    TEST_EQ(0.0, q.a[1]);
    TEST_EQ(0.0, q.a[2]);
    TEST_EQ(0.0, q.a[3]);
    q = make_qdouble("12e+0");
    TEST_EQ(12, q.a[0]);
    TEST_EQ(0.0, q.a[1]);
    TEST_EQ(0.0, q.a[2]);
    TEST_EQ(0.0, q.a[3]);
    q = make_qdouble("12e-0");
    TEST_EQ(12, q.a[0]);
    TEST_EQ(0.0, q.a[1]);
    TEST_EQ(0.0, q.a[2]);
    TEST_EQ(0.0, q.a[3]);
    q = make_qdouble("12e");
    TEST_EQ(12, q.a[0]);
    TEST_EQ(0.0, q.a[1]);
    TEST_EQ(0.0, q.a[2]);
    TEST_EQ(0.0, q.a[3]);
    }
  
  }


void run_all_qdouble_tests()
  {
  test_two_sum_1();
  test_two_sum_2();
  test_two_sum_3();
  test_two_sum_4();
  test_simple_arithmetic();
  test_output_test();  
  test_polynomial();  
  test_machin_formula_pi();   
  test_salamin_brent_quadratic_formula_pi(); 
  test_borwein_quartic_formula_for_pi();  
  test_taylor_series_for_e();   
  test_taylor_series_for_log_2();  
  test_sanity_check_exp();  
  test_sanity_check_sin_cos(); 
  test_make_qdouble();
  }