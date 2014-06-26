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
#include <math.h>
#include <time.h>

#include "bwa.h"
#include "bwamem.h"
#include "bwtindex.h"
#include "kvec.h"
#include "kseq.h"
#include "utils.h"
KSEQ_DECLARE(gzFile)

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

// Modified based on bwa_index in bwtindex.c
int libbwa_index(const char *db, const char *prefix_, libbwa_index_algo algo, int is_64)
{
    char *prefix, *str, *str2, *str3;
    clock_t t;
    int64_t l_pac;

    if (is_64 < 0 || 1 < is_64) return 1;

    prefix = (char *)calloc(strlen(prefix_) + 10, 1);
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

// Same as mem_opt_init in bwamem.c
libbwa_mem_opt *libbwa_mem_opt_init(void)
{
    libbwa_mem_opt *o;
    o = calloc(1, sizeof(mem_opt_t));
    o->flag = 0;
    o->a = 1; o->b = 4;
    o->o_del = o->o_ins = 6;
    o->e_del = o->e_ins = 1;
    o->w = 100;
    o->T = 30;
    o->zdrop = 100;
    o->pen_unpaired = 17;
    o->pen_clip5 = o->pen_clip3 = 5;
    o->min_seed_len = 19;
    o->split_width = 10;
    o->max_occ = 500;
    o->max_chain_gap = 10000;
    o->max_ins = 10000;
    o->mask_level = 0.50;
    o->drop_ratio = 0.50;
    o->XA_drop_ratio = 0.80;
    o->split_factor = 1.5;
    o->chunk_size = 10000000;
    o->n_threads = 1;
    o->max_hits = 5;
    o->max_matesw = 50;
    o->mask_level_redun = 0.95;
    o->min_chain_weight = 0;
    o->max_chain_extend = 1<<30;
    o->mapQ_coef_len = 50; o->mapQ_coef_fac = log(o->mapQ_coef_len);
    bwa_fill_scmat(o->a, o->b, o->mat);
    return o;
}

void convert_mem_opt(const libbwa_mem_opt *src, mem_opt_t *dst)
{
    dst->a = src->a; dst->b = src->b;
    dst->o_del = src->o_del; dst->e_del = src->e_del;
    dst->o_ins = src->o_ins; dst->e_ins = src->e_ins;
    dst->pen_unpaired = src->pen_unpaired;
    dst->pen_clip5 = src->pen_clip5; dst->pen_clip3 = src->pen_clip3;
    dst->w = src->w;
    dst->zdrop = src->zdrop;

    dst->T = src->T;
    dst->flag = src->flag;
    dst->min_seed_len = src->min_seed_len;
    dst->min_chain_weight = src->min_chain_weight;
    dst->max_chain_extend = src->max_chain_extend;
    dst->split_factor = src->split_factor;
    dst->split_width = src->split_width;
    dst->max_occ = src->max_occ;
    dst->max_chain_gap = src->max_chain_gap;
    dst->n_threads = src->n_threads;
    dst->chunk_size = src->chunk_size;
    dst->mask_level = src->mask_level;
    dst->drop_ratio = src->drop_ratio;
    dst->XA_drop_ratio = src->XA_drop_ratio;
    dst->mask_level_redun = src->mask_level_redun;
    dst->mapQ_coef_len = src->mapQ_coef_len;
    dst->mapQ_coef_fac = src->mapQ_coef_fac;
    dst->max_ins = src->max_ins;
    dst->max_matesw = src->max_matesw;
    dst->max_hits = src->max_hits;
    memcpy(dst->mat, src->mat, sizeof(int8_t) * 25);
}

// Modified based on main_mem in fastmap.c
int libbwa_mem(const char *db, const char *reads, const char *mates, const libbwa_mem_opt *opt_)
{
    mem_opt_t *opt;
    int fd, fd2, i, n, copy_comment = 0;
    gzFile fp, fp2 = 0;
    kseq_t *ks, *ks2 = 0;
    bseq1_t *seqs;
    bwaidx_t *idx;
    char *rg_line = 0;
    void *ko = 0, *ko2 = 0;
    int64_t n_processed = 0;
    mem_pestat_t *pes0 = 0;

    opt = mem_opt_init();
    convert_mem_opt(opt_, opt);

    bwa_fill_scmat(opt->a, opt->b, opt->mat);
    if ((idx = bwa_idx_load(db, BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak

    ko = kopen(reads, &fd);
    if (ko == 0) {
        if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, reads);
        return 1;
    }
    fp = gzdopen(fd, "r");
    ks = kseq_init(fp);
    if (mates != NULL) {
        if (opt->flag&MEM_F_PE) {
            if (bwa_verbose >= 2)
                fprintf(stderr, "[W::%s] when '-p' is in use, the second query file will be ignored.\n", __func__);
        } else {
            ko2 = kopen(mates, &fd2);
            if (ko2 == 0) {
                if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, mates);
                return 1;
            }
            fp2 = gzdopen(fd2, "r");
            ks2 = kseq_init(fp2);
            opt->flag |= MEM_F_PE;
        }
    }
    bwa_print_sam_hdr(idx->bns, rg_line);
    while ((seqs = bseq_read(opt->chunk_size * opt->n_threads, &n, ks, ks2)) != 0) {
        int64_t size = 0;
        if ((opt->flag & MEM_F_PE) && (n&1) == 1) {
            if (bwa_verbose >= 2)
                fprintf(stderr, "[W::%s] odd number of reads in the PE mode; last read dropped\n", __func__);
            n = n>>1<<1;
        }
        if (!copy_comment)
            for (i = 0; i < n; ++i) {
                free(seqs[i].comment); seqs[i].comment = 0;
            }
        for (i = 0; i < n; ++i) size += seqs[i].l_seq;
        if (bwa_verbose >= 3)
            fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, n, (long)size);
        mem_process_seqs(opt, idx->bwt, idx->bns, idx->pac, n_processed, n, seqs, pes0);
        n_processed += n;
        for (i = 0; i < n; ++i) {
            err_fputs(seqs[i].sam, stdout);
            free(seqs[i].name); free(seqs[i].comment); free(seqs[i].seq); free(seqs[i].qual); free(seqs[i].sam);
        }
        free(seqs);
    }

    free(opt);
    bwa_idx_destroy(idx);
    kseq_destroy(ks);
    err_gzclose(fp); kclose(ko);
    if (ks2) {
        kseq_destroy(ks2);
        err_gzclose(fp2); kclose(ko2);
    }
    return 0;
}
