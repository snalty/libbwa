libbwa
======

libbwa is a shared library of Burrows-Wheeler Aligner (BWA) forked from original
[repository][bwa].

[![Build Status](https://travis-ci.org/chrovis/libbwa.svg?branch=master)](https://travis-ci.org/chrovis/libbwa)

BWA is a great software in bioinformatics which efficiently maps sequences
against a large reference genome. BWA is, however, implemented as a command-line
tool, and therefore, it is hard-to-use in programs and it lacks reusability. The
purpose of libbwa is providing a shared library of BWA for calling BWA functions
from your programs.

To tracking original BWA changes, libbwa is designed in order that original
codes are not modified as much as possible. In addition, it uses CMake as a
build tool instead of Make for general platforms and prepares unit tests for
ensuring the software safety.

Requirements
------------

- CMake
- zlib

Installation
------------

```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make install
```

Usage
-----

Include `libbwa.h` in your source and compile it with `libbwa` library.
libbwa wraps main BWA commands such as `index`, `mem`, etc.

```c
#include <libbwa.h>

void index_example(void)
{
    char *db = "path/to/reference.fa";
    char *prefix = "path/to/reference.fa";

    // Equivalent to `bwa index` command
    libbwa_index(db, prefix, LIBBWA_INDEX_ALGO_AUTO, 0);
}

void mem_example(void)
{
    char *db = "path/to/reference.fa";
    char *read = "path/to/read.fq";
    char *out = "path/to/out.sam";

    // Create option struct
    libbwa_mem_opt *opt = libbwa_mem_opt_init();

    // Equivalent to `bwa mem` command
    libbwa_mem(db, read, NULL, out, opt);

    // Release option
    libbwa_mem_opt_destroy(opt);
}
```

Unfortunately libbwa reference has not been provided yet.
Please check [libbwa.h][libbwa.h] for more information.

Tests
-----

[CUnit][cunit] is required to run tests.
If you want to run tests, `cmake` with `BUILD_TESTING` option.

```bash
$ cmake -DBUILD_TESTING=ON ..
$ make
$ make test
```

License
-------

Copyright 2014 [Xcoo, Inc.][xcoo]

libbwa is released under [GPLv3][gplv3].
libbwa is based heavily on [BWA][bwa] by [Heng Li][lh3].

[bwa]: https://github.com/lh3/bwa
[libbwa.h]: https://github.com/chrovis/libbwa/blob/master/libbwa.h
[cunit]: http://cunit.sourceforge.net/
[xcoo]: http://www.xcoo.jp/
[gplv3]: http://www.gnu.org/licenses/gpl-3.0.html
[lh3]: https://github.com/lh3
