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
#define LIBBWA_PACKAGE_VERSION "0.7.10-r789"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * libbwa error codes.
 *
 * libbwa function returns a status code, and you can check the success or cause
 * of the failure. This error codes are newly added, and therefore, note that a
 * part of original sources returns other error codes not existing in the enum
 * below.
 */
// NOTE: Should never use 1 or -1 as an entity of a error code. These magic
//       numbers are used in original BWA codes, and they have conflict risk.
typedef enum {
    LIBBWA_E_SUCCESS = 0,

    // General errors
    LIBBWA_E_INVALID_ARGUMENT = 10,
    LIBBWA_E_INVALID_OPTION = 11,
    LIBBWA_E_FILE_ERROR = 12,

    // BWA errors
    LIBBWA_E_INDEX_ERROR = 20,
    LIBBWA_E_UNMATCHED_SAI = 21,
    LIBBWA_E_NOT_IMPLEMENTATION = 22
} libbwa_error_code;

// index
// --------------------

/**
 * Index generation algorithms.
 *
 * @see libbwa_index()
 */
typedef enum {
    LIBBWA_INDEX_ALGO_AUTO = 0,
    LIBBWA_INDEX_ALGO_DIV = 1,
    LIBBWA_INDEX_ALGO_BWTSW = 2,
    LIBBWA_INDEX_ALGO_IS = 3
} libbwa_index_algo;

/**
 * Index database sequences in the FASTA format.
 *
 * Equivalent to `bwa index`.
 *
 * @see libbwa_index_algo
 */
int libbwa_index(const char *db, const char *prefix_, libbwa_index_algo algo, int is_64);

// aln
// --------------------

/**
 * Option structure for aln function.
 *
 * @see libbwa_aln_opt_init()
 * @see libbwa_aln_opt_destroy()
 * @see libbwa_aln()
 */
// Based on gap_opt_t in bwtaln.h
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

/**
 * Returns initialized libbwa_aln_opt.
 *
 * This function dynamically allocates memory. You need to free the memory by
 * libbwa_aln_opt_destroy() after the process.
 *
 * @see libbwa_aln_opt
 * @see libbwa_aln_opt_destroy()
 */
// Based on gap_init_opt in bwtaln.h
libbwa_aln_opt *libbwa_aln_opt_init(void);

/**
 * Destroy libbwa_aln_opt.
 *
 * @see libbwa_aln_opt
 * @see libbwa_aln_opt_init()
 */
void libbwa_aln_opt_destroy(libbwa_aln_opt *opt);

/**
 * Find the SA coordinates of the input reads.
 *
 * Equivalent to `bwa aln`.
 *
 * @see libbwa_aln_opt
 */
int libbwa_aln(const char *db, const char *read, const char *out,
               const libbwa_aln_opt *opt_);

// samse
// --------------------

/**
 * Option structure for samse function.
 *
 * @see libbwa_samse_opt_init()
 * @see libbwa_samse_opt_destroy()
 * @see libbwa_samse()
 */
typedef struct {
    int n_occ;
    char *rg_line;
} libbwa_samse_opt;

/**
 * Returns initialized libbwa_samse_opt.
 *
 * This function dynamically allocates memory. You need to free the memory by
 * libbwa_samse_opt_destroy() after the process.
 *
 * @see libbwa_samse_opt
 * @see libbwa_samse_opt_destroy()
 */
libbwa_samse_opt *libbwa_samse_opt_init(void);

/**
 * Destroy libbwa_samse_opt.
 *
 * @see libbwa_samse_opt
 * @see libbwa_samse_opt_init()
 */
void libbwa_samse_opt_destroy(libbwa_samse_opt *opt);

/**
 * Generate alignments in the SAM format given single-end reads.
 *
 * Equivalent to `bwa samse`.
 *
 * @see libbwa_samse_opt
 */
int libbwa_samse(const char *db, const char *sai, const char *read,
                 const char *out, const libbwa_samse_opt *opt);

// sampe
// --------------------

/**
 * Option structure for sampe function.
 *
 * @see libbwa_sampe_opt_init()
 * @see libbwa_sampe_opt_destroy()
 * @see libbwa_sampe()
 */
// Based on pe_opt_t in bwtaln.h
typedef struct {
    int max_isize, force_isize;
    int max_occ;
    int n_multi, N_multi;
    int type, is_sw, is_preload;
    double ap_prior;
    char *rg_line;
} libbwa_sampe_opt;

/**
 * Returns initialized libbwa_sampe_opt.
 *
 * This function dynamically allocates memory. You need to free the memory by
 * libbwa_sampe_opt_destroy() after the process.
 *
 * @see libbwa_sampe_opt
 * @see libbwa_sampe_opt_destroy()
 */
// Based on bwa_init_pe_opt in bwtaln.h
libbwa_sampe_opt *libbwa_sampe_opt_init(void);

/**
 * Destroy libbwa_sampe_opt.
 *
 * @see libbwa_sampe_opt
 * @see libbwa_sampe_opt_init()
 */
void libbwa_sampe_opt_destroy(libbwa_sampe_opt *opt);

/**
 * Generate alignments in the SAM format given paired-end reads.
 *
 * Equivalent to `bwa sampe`.
 *
 * @see libbwa_sampe_opt
 */
int libbwa_sampe(const char *db, const char *sai1, const char *sai2,
                 const char *read1, const char *read2, const char *out,
                 const libbwa_sampe_opt *opt);

// bwasw
// --------------------

/**
 * Option structure for sw function.
 *
 * @see libbwa_sw_opt_init()
 * @see libbwa_sw_opt_destroy()
 * @see libbwa_sw()
 */
// Based on bsw2opt_t in bwtsw2.h
typedef struct {
    int skip_sw:8, cpy_cmt:8, hard_clip:16;
    int a, b, q, r, t, qr, bw, max_ins, max_chain_gap;
    int z, is, t_seeds, multi_2nd;
    float mask_level, coef;
    int n_threads, chunk_size;
} libbwa_sw_opt;

/**
 * Returns initialized libbwa_sw_opt.
 *
 * This function dynamically allocates memory. You need to free the memory by
 * libbwa_sw_opt_destroy() after the process.
 *
 * @see libbwa_sw_opt
 * @see libbwa_sw_opt_destroy()
 */
// Based on bsw2_init_opt in bwtsw2.h
libbwa_sw_opt *libbwa_sw_opt_init(void);

/**
 * Destroy libbwa_sw_opt.
 *
 * @see libbwa_sw_opt
 * @see libbwa_sw_opt_init()
 */
void libbwa_sw_opt_destroy(libbwa_sw_opt *opt);

/**
 * Align query sequences in the FASTQ file with BWA-SW algorithm.
 *
 * Equivalent to `bwa bwasw`.
 *
 * @see libbwa_sw_opt
 */
int libbwa_sw(const char *db, const char *read, const char *mate,
              const char *out, const libbwa_sw_opt *opt_);

// mem
// --------------------

/**
 * Option structure for mem function.
 *
 * @see libbwa_mem_opt_init()
 * @see libbwa_mem_opt_destroy()
 * @see libbwa_mem()
 */
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

/**
 * Returns initialized libbwa_mem_opt.
 *
 * This function dynamically allocates memory. You need to free the memory by
 * libbwa_mem_opt_destroy() after the process.
 *
 * @see libbwa_mem_opt
 * @see libbwa_mem_opt_destroy()
 */
// Same as mem_opt_init in bwamem.h
libbwa_mem_opt *libbwa_mem_opt_init(void);

/**
 * Destroy libbwa_mem_opt.
 *
 * @see libbwa_mem_opt
 * @see libbwa_mem_opt_init()
 */
void libbwa_mem_opt_destroy(libbwa_mem_opt *opt);

/**
 * Align 70bp-1Mbp query sequences with the BWA-MEM algorithm.
 *
 * Equivalent to `bwa mem`.
 *
 * @see libbwa_mem_opt
 */
int libbwa_mem(const char *db, const char *read, const char *mate,
               const char *out, const libbwa_mem_opt *opt_);

// fastmap
// --------------------

/**
 * Option structure for fastmap function.
 *
 * @see libbwa_fastmap_opt_init()
 * @see libbwa_fastmap_opt_destroy()
 * @see libbwa_fastmap()
 */
typedef struct {
    int print_seq;
    int min_iwidth;
    int min_len;
} libbwa_fastmap_opt;

/**
 * Returns initialized libbwa_fastmap_opt.
 *
 * This function dynamically allocates memory. You need to free the memory by
 * libbwa_fastmap_opt_destroy() after the process.
 *
 * @see libbwa_fastmap_opt
 * @see libbwa_fastmap_opt_destroy()
 */
libbwa_fastmap_opt *libbwa_fastmap_opt_init(void);

/**
 * Destroy libbwa_fastmap_opt.
 *
 * @see libbwa_fastmap_opt
 * @see libbwa_fastmap_opt_init()
 */
void libbwa_fastmap_opt_destroy(libbwa_fastmap_opt *opt);

/**
 * Identifies super-maximal exact matches.
 *
 * Equivalent to `bwa fastmap`.
 *
 * @see libbwa_fastmap_opt
 */
int libbwa_fastmap(const char *db, const char *read, const char *out,
                   const libbwa_fastmap_opt *opt);

// pemerge
// --------------------

// TODO

// fa2pac
// --------------------

/**
 * Converts FASTA to PAC format.
 *
 * Equivalent to `bwa fa2pac`.
 */
int libbwa_fa2pac(const char *db, const char *prefix, int for_only);

// pac2bwt
// --------------------

/**
 * Generates BWT from PAC.
 *
 * Equivalent to `bwa pac2bwt`.
 */
int libbwa_pac2bwt(const char *pac, const char *out, int use_is);

// pac2bwtgen
// --------------------

/**
 * Alternative algorithm for generating BWT.
 *
 * Equivalent to `bwa pac2bwtgen`.
 */
int libbwa_bwtgen(const char *pac, const char *out);

// bwtupdate
// --------------------

/**
 * Updates .bwt to the new format.
 *
 * Equivalent to `bwa bwtupdate`.
 */
int libbwa_bwtupdate(const char *bwt_);

// bwt2sa
// --------------------

/**
 * Generates SA from BWT and Occ.
 *
 * Equivalent to `bwa bwt2sa`.
 */
int libbwa_bwt2sa(const char *bwt_, const char *out, int sa_intv);

#ifdef __cplusplus
}
#endif

#endif
