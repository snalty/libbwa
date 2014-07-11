/*
 * Copyright (C) 2014  Xcoo, Inc.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * This source implements public interfaces for using BWA functions from codes.
 * Some codes included in this source are based on the original BWA codes.
 */

#include "libbwa.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <zlib.h>

#include "bwtindex.h"
#include "bntseq.h"
#include "utils.h"

// Based on bwa_bwtupdate in bwtindex.c
int libbwa_bwtupdate(const char *bwt_)
{
    bwt_t *bwt;
	bwt = bwt_restore_bwt(bwt_);
	bwt_bwtupdate_core(bwt);
	bwt_dump_bwt(bwt_, bwt);
	bwt_destroy(bwt);
	return 0;
}

// Based on bwa_bwt2sa in bwtindex.c
int libbwa_bwt2sa(const char *bwt_, const char *out, int sa_intv)
{
	bwt_t *bwt;
	bwt = bwt_restore_bwt(bwt_);
	bwt_cal_sa(bwt, sa_intv);
	bwt_dump_sa(out, bwt);
	bwt_destroy(bwt);
	return 0;
}

// Modified based on bwa_index in bwtindex.c
int libbwa_index(const char *db, const char *prefix_, libbwa_index_algo algo, int is_64)
{
    char *prefix, *str, *str2, *str3;
    clock_t t;
    int64_t l_pac;

    if (is_64 < 0 || 1 < is_64) return 1;

    prefix = (char*)calloc(strlen(prefix_) + 10, 1);
    strcpy(prefix, prefix_);
    if (is_64) strcat(prefix, ".64");
    str  = (char*)calloc(strlen(prefix) + 10, 1);
    str2 = (char*)calloc(strlen(prefix) + 10, 1);
    str3 = (char*)calloc(strlen(prefix) + 10, 1);
    printf("2");
    { // nucleotide indexing
        gzFile fp = xzopen(db, "r");
        t = clock();
        fprintf(stderr, "[bwa_index] Pack FASTA... ");
        l_pac = bns_fasta2bntseq(fp, prefix, 0);
        fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
        err_gzclose(fp);
    }
    if (algo == LIBBWA_INDEX_ALGO_AUTO)
        algo = l_pac > 50000000 ? LIBBWA_INDEX_ALGO_BWTSW : LIBBWA_INDEX_ALGO_IS; // set the algorithm for generating BWT
    {
        strcpy(str, prefix); strcat(str, ".pac");
        strcpy(str2, prefix); strcat(str2, ".bwt");
        t = clock();
        fprintf(stderr, "[bwa_index] Construct BWT for the packed sequence...\n");
        if (algo == LIBBWA_INDEX_ALGO_BWTSW) bwt_bwtgen(str, str2);
        else if (algo == LIBBWA_INDEX_ALGO_DIV || algo == LIBBWA_INDEX_ALGO_IS) {
            bwt_t *bwt;
            bwt = bwt_pac2bwt(str, algo == LIBBWA_INDEX_ALGO_IS);
            bwt_dump_bwt(str2, bwt);
            bwt_destroy(bwt);
        }
        fprintf(stderr, "[bwa_index] %.2f seconds elapse.\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    {
        bwt_t *bwt;
        strcpy(str, prefix); strcat(str, ".bwt");
        t = clock();
        fprintf(stderr, "[bwa_index] Update BWT... ");
        bwt = bwt_restore_bwt(str);
        bwt_bwtupdate_core(bwt);
        bwt_dump_bwt(str, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    {
        gzFile fp = xzopen(db, "r");
        t = clock();
        fprintf(stderr, "[bwa_index] Pack forward-only FASTA... ");
        l_pac = bns_fasta2bntseq(fp, prefix, 1);
        fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
        err_gzclose(fp);
    }
    {
        bwt_t *bwt;
        strcpy(str, prefix); strcat(str, ".bwt");
        strcpy(str3, prefix); strcat(str3, ".sa");
        t = clock();
        fprintf(stderr, "[bwa_index] Construct SA from BWT and Occ... ");
        bwt = bwt_restore_bwt(str);
        bwt_cal_sa(bwt, 32);
        bwt_dump_sa(str3, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    free(str3); free(str2); free(str);
    return 0;
}
