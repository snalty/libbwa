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
 * This header provides public interfaces for using BWA functions from codes.
 * Some definitions included in this header are based on the original BWA codes.
 */

#ifndef LIBBWA_H
#define LIBBWA_H

#include <stdio.h>
#include <stdint.h>

#include "bntseq.h"

#define LIBBWA_PG_ID "bwa"
#define LIBBWA_PG_PN "libbwa"
#define LIBBWA_PACKAGE_VERSION "0.7.9a-r786"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    LIBBWA_INDEX_ALGO_AUTO = 0,
    LIBBWA_INDEX_ALGO_DIV = 1,
    LIBBWA_INDEX_ALGO_BWTSW = 2,
    LIBBWA_INDEX_ALGO_IS = 3
} libbwa_index_algo;

int libbwa_index(const char *db, const char *prefix_, libbwa_index_algo algo, int is_64);

// Based on pe_opt_t in bwtaln.h
typedef struct {
	int s_mm, s_gapo, s_gape;
	int mode; // bit 24-31 are the barcode length
	int indel_end_skip, max_del_occ, max_entries;
	float fnr;
	int max_diff, max_gapo, max_gape;
	int max_seed_diff, seed_len;
	int n_threads;
	int max_top2;
	int trim_qual;
} libbwa_aln_opt;

// Based on gap_init_opt in bwtaln.h
libbwa_aln_opt *libbwa_aln_opt_init(void);

int libbwa_aln(const char *db, const char *read, const char *out, const libbwa_aln_opt *opt_);

typedef struct {
    int n_occ;
    char *rg_line;
} libbwa_samse_opt;

libbwa_samse_opt *libbwa_samse_opt_init(void);

int libbwa_samse(const char *db, const char *sai, const char *read, const char *out, const libbwa_samse_opt *opt);

// Based on bsw2opt_t in bwtsw2.h
typedef struct {
    int skip_sw:8, cpy_cmt:8, hard_clip:16;
    int a, b, q, r, t, qr, bw, max_ins, max_chain_gap;
    int z, is, t_seeds, multi_2nd;
    float mask_level, coef;
    int n_threads, chunk_size;
} libbwa_sw_opt;

// Based on bsw2_init_opt in bwtsw2.h
libbwa_sw_opt *libbwa_sw_opt_init(void);

int libbwa_sw(const char *db, const char *read, const char *mate, const char *out, const libbwa_sw_opt *opt_);

// Same as mem_opt_t in bwamem.h
typedef struct {
    int a, b;               // match score and mismatch penalty
    int o_del, e_del;
    int o_ins, e_ins;
    int pen_unpaired;       // phred-scaled penalty for unpaired reads
    int pen_clip5,pen_clip3;// clipping penalty. This score is not deducted from the DP score.
    int w;                  // band width
    int zdrop;              // Z-dropoff

    int T;                  // output score threshold; only affecting output
    int flag;               // see MEM_F_* macros
    int min_seed_len;       // minimum seed length
    int min_chain_weight;
    int max_chain_extend;
    float split_factor;     // split into a seed if MEM is longer than min_seed_len*split_factor
    int split_width;        // split into a seed if its occurence is smaller than this value
    int max_occ;            // skip a seed if its occurence is larger than this value
    int max_chain_gap;      // do not chain seed if it is max_chain_gap-bp away from the closest seed
    int n_threads;          // number of threads
    int chunk_size;         // process chunk_size-bp sequences in a batch
    float mask_level;       // regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
    float drop_ratio;       // drop a chain if its seed coverage is below drop_ratio times the seed coverage of a better chain overlapping with the small chain
    float XA_drop_ratio;    // when counting hits for the XA tag, ignore alignments with score < XA_drop_ratio * max_score; only effective for the XA tag
    float mask_level_redun;
    float mapQ_coef_len;
    int mapQ_coef_fac;
    int max_ins;            // when estimating insert size distribution, skip pairs with insert longer than this value
    int max_matesw;         // perform maximally max_matesw rounds of mate-SW for each end
    int max_hits;           // if there are max_hits or fewer, output them all
    int8_t mat[25];         // scoring matrix; mat[0] == 0 if unset
} libbwa_mem_opt;

// Same as mem_opt_init in bwamem.h
libbwa_mem_opt *libbwa_mem_opt_init(void);

int libbwa_mem(const char *db, const char *read, const char *mate, const char *out, const libbwa_mem_opt *opt_);

#ifdef __cplusplus
}
#endif

#endif
