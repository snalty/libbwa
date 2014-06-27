libbwa
======

libbwa is a shared library of Burrows-Wheeler Aligner (BWA) forked from original [repository][bwa].

Requirements
------------

- CMake
- Clang
- zlib
- CUnit (optional)

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

Current libbwa only supports 'index' and 'mem' functions.

```c
# include <libbwa.h>

void index_example(void)
{
    char *db = "path/to/file.fa";
    char *prefix = "path/to/file.fa";

    libbwa_index(db, prefix, LIBBWA_INDEX_ALGO_AUTO, 0);
}

void mem_example(void)
{
    char *db = "path/to/file.fa";
    char *read = "path/to/file.fq";
    char *out = "path/to/out.sam";
    libbwa_mem_opt *opt = libbwa_mem_opt_init();

    libbwa_mem(db, read, NULL, out, opt);
}
```

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
[cunit]: http://cunit.sourceforge.net/
[xcoo]: http://www.xcoo.jp/
[gplv3]: http://www.gnu.org/licenses/gpl-3.0.html
[lh3]: https://github.com/lh3
