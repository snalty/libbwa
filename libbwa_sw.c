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

#include <stdlib.h>

#include "bwtsw2.h"
#include "bwa.h"

// Based on bsw2_init_opt in bwtsw2_aux.c
libbwa_sw_opt *libbwa_sw_opt_init(void)
{
    libbwa_sw_opt *o = (libbwa_sw_opt *)calloc(1, sizeof(libbwa_sw_opt));
    o->a = 1;
    o->b = 3;
    o->q = 5;
    o->r = 2;
    o->t = 30;

    o->bw = 50;

    o->max_ins = 20000;

    o->z = 1;
    o->is = 3;
    o->t_seeds = 5;
    o->hard_clip = 0;
    o->skip_sw = 0;

    o->mask_level = 0.50f;
    o->coef = 5.5f;

    o->qr = o->q + o->r;
    o->n_threads = 1;
    o->chunk_size = 10000000;

    o->max_chain_gap = 10000;

    o->cpy_cmt = 0;

    return o;
}

void convert_sw_opt(const libbwa_sw_opt *src, bsw2opt_t *dst)
{
    dst->a = src->a;
    dst->b = src->b;
    dst->q = src->q;
    dst->r = src->r;
    dst->t = src->t;

    dst->bw = src->bw;

    dst->max_ins = src->max_ins;

    dst->z = src->z;
    dst->is = src->is;
    dst->t_seeds = src->t_seeds;
    dst->hard_clip = src->hard_clip;
    dst->skip_sw = src->skip_sw;

    dst->mask_level = src->mask_level;
    dst->coef = src->coef;

    dst->qr = src->qr;
    dst->n_threads = src->n_threads;
    dst->chunk_size = src->chunk_size;

    dst->max_chain_gap = src->max_chain_gap;

    dst->cpy_cmt = src->cpy_cmt;
}

// Based on bwa_bwtsw2 in bwtsw2_main.c
int libbwa_sw(const char *db, const char *read, const char *mate, const libbwa_sw_opt *opt_)
{
    bsw2opt_t *opt;
    bwaidx_t *idx;

    opt = bsw2_init_opt();
    convert_sw_opt(opt_, opt);

    srand48(11);

    opt->t *= opt->a;
    opt->coef *= opt->a;

    if ((idx = bwa_idx_load(db, BWA_IDX_BWT|BWA_IDX_BNS)) == 0) return 1;
    bsw2_aln(opt, idx->bns, idx->bwt, read, mate != NULL ? mate : 0);
    bwa_idx_destroy(idx);
    free(opt);

    return 0;
}
