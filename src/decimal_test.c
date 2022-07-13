#include <check.h>
#include <math.h>
#include <stdlib.h>

#include "s21_decimal.h"

#define NaN 0.0 / 0.0
#define inf 1.0 / 0.0
#define neg_inf -1.0 / 0.0

START_TEST(s21_add_test) {
  s21_decimal case_full1 = {{123, 321, -2147483648, 0}};
  s21_decimal case_full2 = {{123, 321, -2147483648, 0}};
  s21_decimal result1 = {{0, 0, 0, 0}};
  s21_decimal result2 = {{0, 0, 0, 0}};
  s21_decimal case_max = {{0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0}};
  s21_decimal case_min = {{0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x80000000}};
  int result_full1 = 0;
  result_full1 = s21_add(case_full1, case_full2, &result1);
  ck_assert_int_eq(result_full1, 1);
  case_full1.sign = 1;
  int result_full3 = 0;
  result_full3 = s21_add(case_full1, case_full2, &result1);
  ck_assert_int_eq(result_full3, 0);
  s21_decimal case_full = {{123, 321, 0xFFFFFFFF, 0}};
  result_full1 = s21_add(case_full, case_full, &result2);
  ck_assert_int_eq(result_full1, 1);
  int result_full2 = 0;
  result_full2 = s21_add(case_max, case_min, &result1);
  ck_assert_int_eq(result_full2, 0);
}
END_TEST

START_TEST(s21_sub_test) {
  s21_decimal result1 = {{0, 0, 0, 0}};
  s21_decimal result2 = {{0, 0, 0, 0}};
  s21_decimal case_full1 = {{123, 321, -2147483648, 0}};
  s21_decimal case_full2 = {{123, 321, -2147483648, 0}};
  s21_decimal case1 = {{0, 0, 0, 0}};
  s21_decimal case2 = {{1, 0, 0, 0}};
  s21_decimal case_max = {{0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0}};

  int int_res = -1;
  int res_full1 = 0;
  int res_full2 = 0;
  s21_decimal case_full = {{123, 321, 0xFFFFFFFF, 0}};
  res_full1 = s21_sub(case_full, case_full, &result1);
  ck_assert_int_eq(res_full1, 0);
  setBit(&case_full1, 127, 1);
  res_full2 = s21_sub(case_full1, case_full2, &result2);
  ck_assert_int_eq(res_full2, 2);
  int_res = s21_sub(case1, case_max, &result1);
  ck_assert_int_eq(int_res, 0);
  s21_sub(case1, case2, &result2);
  s21_from_decimal_to_int(result2, &int_res);
  ck_assert_int_eq(int_res, -1);
}
END_TEST

START_TEST(s21_mul_test) {
  s21_decimal result1 = {{0, 0, 0, 0}};
  s21_decimal result2 = {{0, 0, 0, 0}};
  s21_decimal case_full1 = {{123, 321, -2147483648, 0}};
  s21_decimal case_full2 = {{123, 321, -2147483648, 0}};
  s21_decimal case_max = {{0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0}};
  int res_mul = -1;
  s21_decimal case_d = {{0, 0, 0, 0}};
  res_mul = s21_mul(case_d, case_max, &result1);
  ck_assert_int_eq(res_mul, 0);
  s21_decimal max = {{0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0}};
  res_mul = s21_mul(max, max, &result1);
  ck_assert_int_eq(res_mul, 1);
  int res_ful1 = 0;
  int res_ful2 = 0;
  res_ful1 = s21_mul(case_full1, case_full2, &result1);
  ck_assert_int_eq(res_ful1, 1);
  res_ful2 = s21_mul(case_full1, case_full2, &result2);
  setBit(&result2, 127, 1);
  ck_assert_int_eq(res_ful2, 1);
}
END_TEST

START_TEST(s21_div_test) {
  s21_decimal result1 = {{0, 0, 0, 0}};
  s21_decimal result2 = {{0, 0, 0, 0}};
  s21_decimal case_full1 = {{123, 321, -2147483648, 0}};
  s21_decimal case_full2 = {{123, 321, -2147483648, 0}};
  s21_decimal case1 = {{0, 0, 0, 0}};
  s21_decimal case_max = {{0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0}};
  s21_decimal case12 = {{150, 0, 0, 0x80000000}};
  s21_decimal case13 = {{2, 0, 0, 0x80000000}};
  s21_decimal d1 = {{1, 0, 0, 0}};
  d1.exp = 14;
  d1.sign = 1;
  int res = -1;
  int int_res = -1;
  int res_div = -1;
  int_res = s21_div(case1, case_max, &result1);
  ck_assert_int_eq(int_res, 0);
  res_div = s21_div(case_max, case1, &result2);
  ck_assert_int_eq(res_div, 3);
  s21_div(case12, case13, &result2);
  int res_ful1 = 0;
  int res_ful2 = 0;
  res_ful1 = s21_div(case_full1, case_full2, &result1);
  ck_assert_int_eq(res_ful1, 0);
  s21_decimal one = {{1, 0, 0, 0}};
  res_ful2 = s21_div(one, case_max, &result2);
  ck_assert_int_eq(res_ful2, 2);
  res_ful1 = s21_mod(one, one, &result1);
  ck_assert_int_eq(res_ful1, 0);
  res = s21_div(case_max, d1, &result1);
  ck_assert_int_eq(res, 2);
}
END_TEST

START_TEST(s21_mod_test) {
  {
    s21_decimal src1, src2, res_mod;
    int a = 3;
    int b = 2;
    int res_origin = a % b;
    int res = 0;
    s21_from_int_to_decimal(a, &src1);
    s21_from_int_to_decimal(b, &src2);
    s21_mod(src1, src2, &res_mod);
    s21_from_decimal_to_int(res_mod, &res);
    ck_assert_int_eq(res_origin, res);
  }
  {
    s21_decimal src1, src2, res_mod;
    float a = 1.2;
    int b = 3;
    float res_origin = fmod(a, b);
    float res = 0;
    s21_from_float_to_decimal(a, &src1);
    s21_from_int_to_decimal(b, &src2);
    s21_mod(src1, src2, &res_mod);
    s21_from_decimal_to_float(res_mod, &res);
    ck_assert_float_eq(res_origin, res);
  }
  {
    s21_decimal src1, src2, res_mod;
    initDecimal(&src1);
    initDecimal(&src2);
    initDecimal(&res_mod);
    double a = 98765;
    double b = -1234;
    int res_origin = (int)fmod(a, b);
    int res = 0;
    s21_from_int_to_decimal(a, &src1);
    s21_from_int_to_decimal(b, &src2);
    s21_mod(src1, src2, &res_mod);
    s21_from_decimal_to_int(res_mod, &res);
    ck_assert_int_eq(res_origin, res);
  }
  {
    s21_decimal src1, src2, res_mod;
    int a = 30198;
    int b = 20210;
    int res_origin = a % b;
    int res = 0;
    s21_from_int_to_decimal(a, &src1);
    s21_from_int_to_decimal(b, &src2);
    s21_mod(src1, src2, &res_mod);
    s21_from_decimal_to_int(res_mod, &res);
    ck_assert_int_eq(res_origin, res);
  }
  {
    s21_decimal src1, src2, res_mod;
    int a = -98765;
    int b = -1234;
    int res_origin = a % b;
    int res = 0;
    s21_from_int_to_decimal(a, &src1);
    s21_from_int_to_decimal(b, &src2);
    s21_mod(src1, src2, &res_mod);
    s21_from_decimal_to_int(res_mod, &res);
    ck_assert_int_eq(res_origin, res);
  }
  {
    s21_decimal src1, src2, res_mod;
    int a = 98765;
    int b = 127234;
    int res_origin = a % b;
    int res = 0;
    s21_from_int_to_decimal(a, &src1);
    s21_from_int_to_decimal(b, &src2);
    s21_mod(src1, src2, &res_mod);
    s21_from_decimal_to_int(res_mod, &res);
    ck_assert_int_eq(res_origin, res);
  }
  {
    s21_decimal src1, src2, res_mod;
    int a = 342576;
    int b = 1542134;
    int res_origin = a % b;
    int res = 0;
    s21_from_int_to_decimal(a, &src1);
    s21_from_int_to_decimal(b, &src2);
    s21_mod(src1, src2, &res_mod);
    s21_from_decimal_to_int(res_mod, &res);
    ck_assert_int_eq(res_origin, res);
  }
  {
    s21_decimal src1, src2, res_mod;
    float a = 1.2;
    float b = 0.3;
    float res_origin = fmod(a, b);
    float res = 0.0;
    s21_from_float_to_decimal(a, &src1);
    s21_from_float_to_decimal(b, &src2);
    s21_mod(src1, src2, &res_mod);
    s21_from_decimal_to_float(res_mod, &res);
    ck_assert_float_eq(res_origin, res);
  }
  {
    s21_decimal src1, src2, res_mod;
    float a = 18;
    float b = -6;
    float res_origin = fmod(a, b);
    float res = 0.0;
    s21_from_float_to_decimal(a, &src1);
    s21_from_float_to_decimal(b, &src2);
    s21_mod(src1, src2, &res_mod);
    s21_from_decimal_to_float(res_mod, &res);
    ck_assert_float_eq(res_origin, res);
  }
  {
    s21_decimal src1 = {{0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0}};
    s21_decimal src2 = {{15, 0, 0, 0}};
    s21_decimal res_mod;
    src2.exp = 3;
    ck_assert_int_eq(s21_mod(src1, src2, &res_mod), 0);
  }
}
END_TEST

START_TEST(s21_is_less_test) {
  s21_decimal test1 = {{0, 0, 0, 0}};
  s21_decimal test2 = {{1, 0, 0, 0}};
  s21_decimal test3 = {{4294967295, 4294967295, 4294967295, 2147483648}};
  s21_decimal test4 = {{4294967295, 4294967295, 4294967295, 0}};

  ck_assert_int_eq(s21_is_less(test1, test2), 1);
  ck_assert_int_eq(s21_is_less(test2, test2), 0);
  ck_assert_int_eq(s21_is_less(test3, test4), 1);
  ck_assert_int_eq(s21_is_less(test4, test3), 0);
  ck_assert_int_eq(s21_is_less(test4, test1), 0);

  s21_decimal test8 = {{0, 0, 1, 0}};
  s21_decimal test9 = {{0, 0, 1, 0}};
  setBit(&test8, 0, 1);
  setBit(&test9, 2, 1);
  ck_assert_int_eq(s21_is_less(test8, test9), 1);
}
END_TEST

START_TEST(s21_is_less_or_equal_test) {
  s21_decimal test1 = {{7, 0, 0, 0}};
  s21_decimal test2 = {{7, 0, 0, 0}};
  s21_decimal test3 = {{4294967295, 4294967295, 4294967295, 0}};
  s21_decimal test4 = {{4294967295, 4294967295, 4294967295, 0}};

  ck_assert_int_eq(s21_is_less_or_equal(test1, test2), 1);
  setBit(&test1, 127, 1);
  ck_assert_int_eq(s21_is_less_or_equal(test1, test2), 1);
  ck_assert_int_eq(s21_is_less_or_equal(test2, test2), 1);
  ck_assert_int_eq(s21_is_less_or_equal(test3, test4), 1);
  setBit(&test3, 127, 1);
  ck_assert_int_eq(s21_is_less_or_equal(test3, test4), 1);
  ck_assert_int_eq(s21_is_less_or_equal(test4, test3), 0);
  ck_assert_int_eq(s21_is_less_or_equal(test4, test1), 0);
  setBit(&test4, 127, 1);
  ck_assert_int_eq(s21_is_less_or_equal(test4, test1), 1);

  s21_decimal test8 = {{0, 0, 1, 0}};
  s21_decimal test9 = {{0, 0, 1, 0}};
  ck_assert_int_eq(s21_is_less_or_equal(test8, test9), 1);
  setBit(&test8, 0, 1);
  setBit(&test9, 0, 1);
  ck_assert_int_eq(s21_is_less_or_equal(test8, test9), 1);
}
END_TEST

START_TEST(s21_is_greater_test) {
  s21_decimal test1 = {{777, 0, 0, 0}};
  s21_decimal test2 = {{777, 0, 0, 0}};
  s21_decimal test3 = {{4294967295, 4294967295, 4294967295, 0}};
  s21_decimal test4 = {{4294967295, 4294967295, 4294967295, 0}};

  ck_assert_int_eq(s21_is_greater(test1, test2), 0);
  setBit(&test3, 127, 1);
  ck_assert_int_eq(s21_is_greater(test1, test2), 0);
  ck_assert_int_eq(s21_is_greater(test2, test2), 0);
  ck_assert_int_eq(s21_is_greater(test3, test4), 0);
  setBit(&test3, 127, 1);
  ck_assert_int_eq(s21_is_greater(test3, test4), 0);
  ck_assert_int_eq(s21_is_greater(test4, test3), 1);
  ck_assert_int_eq(s21_is_greater(test1, test4), 0);

  s21_decimal test8 = {{0, 0, 1, 0}};
  s21_decimal test9 = {{0, 0, 1, 0}};
  setBit(&test8, 0, 1);
  setBit(&test9, 0, 1);
  ck_assert_int_eq(s21_is_greater(test8, test9), 0);
}
END_TEST

START_TEST(s21_is_greater_or_equal_test) {
  s21_decimal test1 = {{777, 777, 777, 0}};
  s21_decimal test2 = {{777, 777, 777, 0}};
  s21_decimal test3 = {{4294967295, 4294967295, 4294967295, 0}};
  s21_decimal test4 = {{4294967295, 4294967295, 4294967295, 0}};

  ck_assert_int_eq(s21_is_greater_or_equal(test1, test2), 1);
  setBit(&test2, 0, 1);
  ck_assert_int_eq(s21_is_greater_or_equal(test1, test2), 1);
  ck_assert_int_eq(s21_is_greater_or_equal(test2, test2), 1);
  ck_assert_int_eq(s21_is_greater_or_equal(test3, test4), 1);
  setBit(&test4, 127, 1);
  ck_assert_int_eq(s21_is_greater_or_equal(test3, test4), 1);
  ck_assert_int_eq(s21_is_greater_or_equal(test4, test3), 0);
  ck_assert_int_eq(s21_is_greater_or_equal(test1, test3), 0);
  setBit(&test4, 127, 1);
  ck_assert_int_eq(s21_is_greater_or_equal(test1, test3), 0);

  s21_decimal test8 = {{0, 0, 1, 0}};
  s21_decimal test9 = {{0, 0, 1, 0}};
  ck_assert_int_eq(s21_is_greater_or_equal(test8, test9), 1);
  setBit(&test8, 0, 1);
  setBit(&test9, 0, 1);
  ck_assert_int_eq(s21_is_greater_or_equal(test8, test9), 1);
}
END_TEST

START_TEST(s21_is_equal_test) {
  s21_decimal test1 = {{777, 777, 777, 0}};
  s21_decimal test2 = {{777, 777, 777, 0}};
  s21_decimal test3 = {{4294967295, 4294967295, 4294967295, 0}};
  s21_decimal test4 = {{4294967295, 4294967295, 4294967295, 0}};

  ck_assert_int_eq(s21_is_equal(test1, test2), 1);
  setBit(&test1, 0, 1);
  ck_assert_int_eq(s21_is_equal(test1, test2), 1);
  ck_assert_int_eq(s21_is_equal(test3, test3), 1);
  ck_assert_int_eq(s21_is_equal(test3, test4), 1);
  setBit(&test4, 127, 1);
  ck_assert_int_eq(s21_is_equal(test3, test4), 0);
}
END_TEST

START_TEST(s21_is_not_equal_test) {
  s21_decimal test1 = {{777, 777, 777, 0}};
  s21_decimal test2 = {{777, 777, 777, 0}};
  s21_decimal test3 = {{4294967295, 4294967295, 4294967295, 0}};
  s21_decimal test4 = {{4294967295, 4294967295, 4294967295, 0}};

  ck_assert_int_eq(s21_is_not_equal(test1, test2), 0);
  setBit(&test4, 0, 1);
  ck_assert_int_eq(s21_is_not_equal(test1, test2), 0);
  ck_assert_int_eq(s21_is_not_equal(test3, test3), 0);
  ck_assert_int_eq(s21_is_not_equal(test3, test4), 0);
  setBit(&test4, 127, 1);
  ck_assert_int_eq(s21_is_not_equal(test3, test4), 1);
}
END_TEST

START_TEST(s21_from_int_to_decimal_test) {
  {
    int test_int = 0;
    s21_decimal *decimal = NULL;
    int status = s21_from_int_to_decimal(test_int, decimal);
    ck_assert_int_eq(status, 1);
  }
  {
    int test_int = 0;
    s21_decimal *decimal = malloc(sizeof(s21_decimal));
    int status = s21_from_int_to_decimal(test_int, decimal);
    ck_assert_int_eq(status, 0);
    ck_assert_int_eq(test_int, decimal->bits[0]);
    ck_assert_int_eq(0, decimal->bits[1]);
    ck_assert_int_eq(0, decimal->bits[2]);
    ck_assert_int_eq(0, getBit(*decimal, 127));
    ck_assert_int_eq(0, getBit(*decimal, 127));
    free(decimal);
  }
  {
    int test_int = INT32_MAX;
    s21_decimal *decimal = malloc(sizeof(s21_decimal));
    int status = s21_from_int_to_decimal(test_int, decimal);
    ck_assert_int_eq(status, 0);
    ck_assert_int_eq(test_int, decimal->bits[0]);
    ck_assert_int_eq(0, decimal->bits[1]);
    ck_assert_int_eq(0, decimal->bits[2]);
    ck_assert_int_eq(0, getBit(*decimal, 127));
    ck_assert_int_eq(0, getBit(*decimal, 127));
    free(decimal);
  }
  {
    int test_int = INT32_MIN + 1;
    s21_decimal *decimal = malloc(sizeof(s21_decimal));
    int status = s21_from_int_to_decimal(-test_int, decimal);
    ck_assert_int_eq(status, 0);
    ck_assert_int_eq(-test_int, decimal->bits[0]);
    ck_assert_int_eq(0, decimal->bits[1]);
    ck_assert_int_eq(0, decimal->bits[2]);
    ck_assert_int_eq(0, getBit(*decimal, 127));
    ck_assert_int_eq(0, getBit(*decimal, 127));
    free(decimal);
  }
  {
    int test_int = -126123;
    s21_decimal *decimal = malloc(sizeof(s21_decimal));
    int status = s21_from_int_to_decimal(-test_int, decimal);
    ck_assert_int_eq(status, 0);
    ck_assert_int_eq(-test_int, decimal->bits[0]);
    ck_assert_int_eq(0, decimal->bits[1]);
    ck_assert_int_eq(0, decimal->bits[2]);
    ck_assert_int_eq(0, getBit(*decimal, 127));
    ck_assert_int_eq(0, getBit(*decimal, 127));
    free(decimal);
  }
  {
    int test_int = 3752;
    s21_decimal *decimal = malloc(sizeof(s21_decimal));
    int status = s21_from_int_to_decimal(test_int, decimal);
    ck_assert_int_eq(status, 0);
    ck_assert_int_eq(test_int, decimal->bits[0]);
    ck_assert_int_eq(0, decimal->bits[1]);
    ck_assert_int_eq(0, decimal->bits[2]);
    ck_assert_int_eq(0, getBit(*decimal, 127));
    ck_assert_int_eq(0, getBit(*decimal, 127));
    free(decimal);
  }
}
END_TEST

START_TEST(s21_from_float_to_decimal_test) {
  float zero_f = 0.0, zero_fl;
  {
    s21_decimal dec;
    int r1 = s21_from_float_to_decimal(zero_f, &dec);
    int r2 = s21_from_decimal_to_float(dec, &zero_fl);
    ck_assert_float_eq_tol(zero_f, zero_fl, 1e-6);
    ck_assert_int_eq(r1, 0);
    ck_assert_int_eq(r2, 0);
  }
  {
    s21_decimal dec;
    int r1 = s21_from_float_to_decimal(1e+29, &dec);
    ck_assert_int_eq(r1, 1);
  }
  {
    s21_decimal dec;
    int r1 = s21_from_float_to_decimal(1e-29, &dec);
    int r2 = s21_from_decimal_to_float(dec, &zero_f);
    ck_assert_int_eq(r1, 1);
    ck_assert_int_eq(r2, 0);
    ck_assert_float_eq(zero_f, 0.0);
  }
  {
    s21_decimal dec;
    int r1 = s21_from_float_to_decimal(NaN, &dec);
    ck_assert_int_eq(r1, 1);
  }
  {
    s21_decimal dec;
    int r1 = s21_from_float_to_decimal(inf, &dec);
    ck_assert_int_eq(r1, 1);
  }
  {
    s21_decimal dec;
    int r1 = s21_from_float_to_decimal(neg_inf, &dec);
    ck_assert_int_eq(r1, 1);
  }
  {
    s21_decimal dec;
    int r1 = s21_from_float_to_decimal(neg_inf, &dec);
    ck_assert_int_eq(r1, 1);
  }
}
END_TEST

START_TEST(s21_from_decimal_to_int_test) {
  s21_decimal dec, big_dec = {{23536, 37375, 5647, 0}};
  int zero_i = 0, n1 = -1233456789, n2 = 98765321;
  int zero_res, res1, res2;
  int r1 = s21_from_int_to_decimal(zero_i, &dec);
  int r2 = s21_from_decimal_to_int(dec, &zero_res);
  ck_assert_int_eq(zero_i, zero_res);
  ck_assert_int_eq(r1, 0);
  ck_assert_int_eq(r2, 0);
  r1 = s21_from_int_to_decimal(n1, &dec);
  r2 = s21_from_decimal_to_int(dec, &res1);
  ck_assert_int_eq(n1, res1);
  ck_assert_int_eq(r1, 0);
  ck_assert_int_eq(r2, 0);
  r1 = s21_from_int_to_decimal(n2, &dec);
  r2 = s21_from_decimal_to_int(dec, &res2);
  ck_assert_int_eq(n2, res2);
  ck_assert_int_eq(r1, 0);
  ck_assert_int_eq(r2, 0);
  r1 = s21_from_int_to_decimal(n2, &dec);
  r2 = s21_from_decimal_to_int(dec, &res2);
  ck_assert_int_eq(n2, res2);
  ck_assert_int_eq(r1, 0);
  ck_assert_int_eq(r2, 0);
  r1 = s21_from_decimal_to_int(big_dec, &res1);
  ck_assert_int_eq(r1, 1);
  big_dec.exp = 15;
  setBit(&big_dec, 127, 1);
  r2 = s21_from_decimal_to_int(big_dec, &res2);
  ck_assert_int_eq(r2, 0);
}
END_TEST

START_TEST(s21_from_decimal_to_float_test) {
  s21_decimal dec;
  float zero_f = 0.0, nmb1 = -654321.123456, nmb2 = 1.678e-10, zero_fl, res1,
        res2;
  int r1 = s21_from_float_to_decimal(nmb1, &dec);
  int r2 = s21_from_decimal_to_float(dec, &res1);
  ck_assert_float_eq_tol(nmb1, res1, 1e-6);
  ck_assert_int_eq(r1, 0);
  ck_assert_int_eq(r2, 0);
  r1 = s21_from_float_to_decimal(nmb2, &dec);
  r2 = s21_from_decimal_to_float(dec, &res2);
  ck_assert_float_eq_tol(nmb2, res2, 1e-6);
  ck_assert_int_eq(r1, 0);
  ck_assert_int_eq(r2, 0);
  r1 = s21_from_float_to_decimal(zero_f, &dec);
  r2 = s21_from_decimal_to_float(dec, &zero_fl);
  ck_assert_float_eq_tol(zero_f, zero_fl, 1e-6);
  ck_assert_int_eq(r1, 0);
  ck_assert_int_eq(r2, 0);
  r1 = s21_from_float_to_decimal(1.36e-20, &dec);
  r2 = s21_from_decimal_to_float(dec, &res1);
  ck_assert_float_eq_tol(1.36e-20, res1, 1e-6);
  ck_assert_int_eq(r1, 0);
  ck_assert_int_eq(r2, 0);
  r1 = s21_from_float_to_decimal(1.36e+20, &dec);
  r2 = s21_from_decimal_to_float(dec, &res1);
  ck_assert_int_eq(r1, 0);
  ck_assert_int_eq(r2, 0);
  r1 = s21_from_float_to_decimal(1e+29, &dec);
  ck_assert_int_eq(r1, 1);
  r1 = s21_from_float_to_decimal(1e-29, &dec);
  r2 = s21_from_decimal_to_float(dec, &zero_f);
  ck_assert_int_eq(r1, 1);
  ck_assert_int_eq(r2, 0);
  ck_assert_float_eq(zero_f, 0.0);
  r1 = s21_from_float_to_decimal(NaN, &dec);
  ck_assert_int_eq(r1, 1);
  r1 = s21_from_float_to_decimal(inf, &dec);
  ck_assert_int_eq(r1, 1);
  r1 = s21_from_float_to_decimal(neg_inf, &dec);
  ck_assert_int_eq(r1, 1);
  r1 = s21_from_float_to_decimal(neg_inf, &dec);
  ck_assert_int_eq(r1, 1);
}
END_TEST

START_TEST(s21_floor_test) {
  float a = 0.4;
  s21_decimal fl_d;
  s21_from_float_to_decimal(a, &fl_d);
  s21_floor(fl_d, &fl_d);
  float result;
  s21_from_decimal_to_float(fl_d, &result);
  ck_assert_double_eq(result, 0.0);

  a = -0.4;
  s21_from_float_to_decimal(a, &fl_d);
  s21_floor(fl_d, &fl_d);
  s21_from_decimal_to_float(fl_d, &result);
  ck_assert_double_eq(result, -1);
  a = 0;
  s21_from_float_to_decimal(a, &fl_d);
  s21_floor(fl_d, &fl_d);
  s21_from_decimal_to_float(fl_d, &result);
  ck_assert_double_eq(result, 0);
}
END_TEST

START_TEST(s21_round_test) {
  float a = 0.4;
  s21_decimal fl_d;
  s21_from_float_to_decimal(a, &fl_d);
  s21_round(fl_d, &fl_d);
  float result;
  s21_from_decimal_to_float(fl_d, &result);
  ck_assert_double_eq(result, 0);
  a = -3.5;
  s21_from_float_to_decimal(a, &fl_d);
  s21_round(fl_d, &fl_d);
  s21_from_decimal_to_float(fl_d, &result);
  ck_assert_double_eq(result, -4);
  s21_decimal d1 = {{12345}};
  d1.exp = 1;
  s21_round(d1, &d1);
  ck_assert_double_eq(d1.bits[0], 1235);
}
END_TEST

START_TEST(s21_truncate_test) {
  float a = 2.6;
  s21_decimal fl_d;
  s21_from_float_to_decimal(a, &fl_d);
  s21_truncate(fl_d, &fl_d);
  float result;
  s21_from_decimal_to_float(fl_d, &result);
  ck_assert_double_eq(result, 2);
  a = -4.5;
  s21_from_float_to_decimal(a, &fl_d);
  s21_truncate(fl_d, &fl_d);
  s21_from_decimal_to_float(fl_d, &result);
  ck_assert_double_eq(result, -4);
}
END_TEST

START_TEST(s21_negate_test) {
  s21_decimal n1 = {{0, 0, 0, 0}};
  s21_decimal n2 = {{0, 0, 0, 0}};
  s21_decimal res = {{0, 0, 0, 0}};
  setBit(&n1, 127, 1);
  s21_negate(n2, &res);
  ck_assert_int_eq(getBit(n1, 127), getBit(res, 127));
}
END_TEST

int main() {
  Suite *suite = suite_create("S21_Debuger");
  SRunner *srunner = srunner_create(suite);

  TCase *s21_add_case = tcase_create("s21_add_case");
  suite_add_tcase(suite, s21_add_case);
  tcase_add_test(s21_add_case, s21_add_test);

  TCase *s21_sub_case = tcase_create("s21_sub_case");
  suite_add_tcase(suite, s21_sub_case);
  tcase_add_test(s21_sub_case, s21_sub_test);

  TCase *s21_mul_case = tcase_create("s21_mul_case");
  suite_add_tcase(suite, s21_mul_case);
  tcase_add_test(s21_mul_case, s21_mul_test);

  TCase *s21_div_case = tcase_create("s21_div_case");
  suite_add_tcase(suite, s21_div_case);
  tcase_add_test(s21_div_case, s21_div_test);

  TCase *s21_mod_case = tcase_create("s21_mod_case");
  suite_add_tcase(suite, s21_mod_case);
  tcase_add_test(s21_mod_case, s21_mod_test);

  TCase *s21_is_less_case = tcase_create("s21_is_less_case");
  suite_add_tcase(suite, s21_is_less_case);
  tcase_add_test(s21_is_less_case, s21_is_less_test);

  TCase *s21_is_less_or_equal_case = tcase_create("s21_is_less_or_equal_case");
  suite_add_tcase(suite, s21_is_less_or_equal_case);
  tcase_add_test(s21_is_less_or_equal_case, s21_is_less_or_equal_test);

  TCase *s21_is_greater_case = tcase_create("s21_is_greater_case");
  suite_add_tcase(suite, s21_is_greater_case);
  tcase_add_test(s21_is_greater_case, s21_is_greater_test);

  TCase *s21_is_greater_or_equal_case =
      tcase_create("s21_is_greater_or_equal_case");
  suite_add_tcase(suite, s21_is_greater_or_equal_case);
  tcase_add_test(s21_is_greater_or_equal_case, s21_is_greater_or_equal_test);

  TCase *s21_is_equal_case = tcase_create("s21_is_equal_case");
  suite_add_tcase(suite, s21_is_equal_case);
  tcase_add_test(s21_is_equal_case, s21_is_equal_test);

  TCase *s21_is_not_equal_case = tcase_create("s21_is_not_equal_case");
  suite_add_tcase(suite, s21_is_not_equal_case);
  tcase_add_test(s21_is_not_equal_case, s21_is_not_equal_test);

  TCase *s21_from_int_to_decimal_case =
      tcase_create("s21_from_int_to_decimal_case");
  suite_add_tcase(suite, s21_from_int_to_decimal_case);
  tcase_add_test(s21_from_int_to_decimal_case, s21_from_int_to_decimal_test);

  TCase *s21_from_float_to_decimal_case =
      tcase_create("s21_from_float_to_decimal_case");
  suite_add_tcase(suite, s21_from_float_to_decimal_case);
  tcase_add_test(s21_from_float_to_decimal_case,
                 s21_from_float_to_decimal_test);

  TCase *s21_from_decimal_to_int_case =
      tcase_create("s21_from_decimal_to_int_case");
  suite_add_tcase(suite, s21_from_decimal_to_int_case);
  tcase_add_test(s21_from_decimal_to_int_case, s21_from_decimal_to_int_test);

  TCase *s21_from_decimal_to_float_case =
      tcase_create("s21_from_decimal_to_float_case");
  suite_add_tcase(suite, s21_from_decimal_to_float_case);
  tcase_add_test(s21_from_decimal_to_float_case,
                 s21_from_decimal_to_float_test);

  TCase *s21_floor_case = tcase_create("s21_floor_case");
  suite_add_tcase(suite, s21_floor_case);
  tcase_add_test(s21_floor_case, s21_floor_test);

  TCase *s21_round_case = tcase_create("s21_round_case");
  suite_add_tcase(suite, s21_round_case);
  tcase_add_test(s21_round_case, s21_round_test);

  TCase *s21_truncate_case = tcase_create("s21_truncate_case");
  suite_add_tcase(suite, s21_truncate_case);
  tcase_add_test(s21_truncate_case, s21_truncate_test);

  TCase *s21_negate_case = tcase_create("s21_negate_case");
  suite_add_tcase(suite, s21_negate_case);
  tcase_add_test(s21_negate_case, s21_negate_test);

  srunner_run_all(srunner, CK_VERBOSE);
  int number_failed = srunner_ntests_failed(srunner);
  srunner_free(srunner);
  return (number_failed == 0) ? 0 : 1;
}
