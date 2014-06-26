/*
 * Copyright 2014 Xcoo, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *
 * This header is added for using two functions, bwa_seq_len and bwt_pac2bwt, in
 * libbwa.
 */

#ifndef BWTINDEX_H
#define BWTINDEX_H

#include "bwt.h"

int64_t bwa_seq_len(const char *fn_pac);
bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is);

#endif
