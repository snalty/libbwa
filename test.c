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
#define TEST_PAC "../test-resources/test.fa.pac"
#define TEST_BWT "../test-resources/test.fa.bwt"

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

void libbwa_sampe_test(void)
{
    char *db = TEST_DB;
    char *sai1 = TEST_SAI, *sai2 = TEST_SAI;
    char *read1 = TEST_READ, *read2 = TEST_READ;
    char *out = "/tmp/libbwa_sampe_test.sam";
    libbwa_sampe_opt *opt = libbwa_sampe_opt_init();
    CU_ASSERT(libbwa_sampe(db, sai1, sai2, read1, read2, out, opt) == 0);
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

void libbwa_fa2pac_test(void)
{
    char *db = TEST_DB;
    char *prefix = "/tmp/libbwa_fa2pac_test.fa";
    int for_only = 0;
    CU_ASSERT(libbwa_fa2pac(db, prefix, for_only) == 0);
}

void libbwa_pac2bwt_test(void)
{
    char *pac = TEST_PAC;
    char *out = "/tmp/libbwa_pac2bwt_test.bwt";
    int use_is = 1;
    CU_ASSERT(libbwa_pac2bwt(pac, out, use_is) == 0);
}

void libbwa_bwtgen_test(void)
{
    char *pac = TEST_PAC;
    char *out = "/tmp/libbwa_bwtgen_test.bwt";
    CU_ASSERT(libbwa_bwtgen(pac, out) == 0);
}

void libbwa_bwtupdate_test(void)
{
    char *bwt = TEST_BWT;
    CU_ASSERT(libbwa_bwtupdate(bwt) == 0);
}

int main(int argc, char *argv[])
{
    unsigned int num_failures;

    CU_TestInfo tests[] = {
        {"index test", libbwa_index_test},
        {"aln test", libbwa_aln_test},
        {"samse test", libbwa_samse_test},
        {"sampe test", libbwa_sampe_test},
        {"sw test", libbwa_sw_test},
        {"mem test", libbwa_mem_test},
        {"fa2pac test", libbwa_fa2pac_test},
        {"pac2bwt test", libbwa_pac2bwt_test},
        {"bwtgen test", libbwa_bwtgen_test},
        // {"bwtupdate test", libbwa_bwtupdate_test},
        CU_TEST_INFO_NULL
    };

    CU_SuiteInfo suites[] = {
        {"Base_Suite", NULL, NULL, tests},
        CU_SUITE_INFO_NULL
    };

   if (CUE_SUCCESS != CU_initialize_registry())
      return CU_get_error();

   if (CUE_SUCCESS != CU_register_suites(suites)) {
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
