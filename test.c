/*
 *  Copyright (C) 2014  Xcoo, Inc.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <CUnit/Basic.h>
#include "libbwa.h"

#define TEST_DB "../test-resources/test.fa"
#define TEST_READ "../test-resources/test.fq"
#define TEST_SAI "../test-resources/test.sai"

char *bwa_pg;

void libbwa_index_test(void)
{
    char *db = TEST_DB;
    char *prefix = "/tmp/test.fa";
    libbwa_index_algo algo = LIBBWA_INDEX_ALGO_AUTO;
    int is_64 = 0;
    CU_ASSERT(libbwa_index(db, prefix, algo, is_64) == 0);
}

void libbwa_aln_test(void)
{
    char *db = TEST_DB;
    char *read = TEST_READ;
    char *out = "/tmp/libbwa_aln_test.sai";
    libbwa_aln_opt *opt = libbwa_aln_opt_init();
    CU_ASSERT(libbwa_aln(db, read, out, opt) == 0);
}

void libbwa_samse_test(void)
{
    char *db = TEST_DB;
    char *sai = TEST_SAI;
    char *read = TEST_READ;
    char *out = "/tmp/libbwa_samse_test.sam";
    libbwa_samse_opt *opt = libbwa_samse_opt_init();
    CU_ASSERT(libbwa_samse(db, sai, read, out, opt) == 0);
}

void libbwa_sw_test(void)
{
    char *db = TEST_DB;
    char *read = TEST_READ;
    char *out = "/tmp/libbwa_sw_test.sam";
    libbwa_sw_opt *opt = libbwa_sw_opt_init();
    CU_ASSERT(libbwa_sw(db, read, NULL, out, opt) == 0);
}

void libbwa_mem_test(void)
{
    char *db = TEST_DB;
    char *read = TEST_READ;
    char *out = "/tmp/libbwa_mem_test.sam";
    libbwa_mem_opt *opt = libbwa_mem_opt_init();
    CU_ASSERT(libbwa_mem(db, read, NULL, out, opt) == 0);
}

int main(int argc, char *argv[])
{
    CU_pSuite pSuite = NULL;
    unsigned int num_failures;

   if (CUE_SUCCESS != CU_initialize_registry())
      return CU_get_error();

   pSuite = CU_add_suite("Base_Suite", NULL, NULL);
   if (NULL == pSuite) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   if (NULL == CU_add_test(pSuite, "index test", libbwa_index_test)) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   if (NULL == CU_add_test(pSuite, "aln test", libbwa_aln_test)) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   if (NULL == CU_add_test(pSuite, "samse test", libbwa_samse_test)) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   if (NULL == CU_add_test(pSuite, "sw test", libbwa_sw_test)) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   if (NULL == CU_add_test(pSuite, "mem test", libbwa_mem_test)) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   CU_basic_set_mode(CU_BRM_VERBOSE);
   CU_basic_run_tests();
   num_failures = CU_get_number_of_failures();
   CU_cleanup_registry();

   if (num_failures > 0) return 1;
   return CU_get_error();
}
