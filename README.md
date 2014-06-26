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

Current libbwa supports only 'index' and 'mem' funtions.

TODO: usage

Tests
-----

If you want to run tests, `cmake` with `BUILD_TESTING` option.

```bash
$ cmake -DBUILD_TESTING=ON ..
$ make
$ make test
```

License
-------

libbwa is released under [GPLv3][gplv3].
libbwa is based heavily on [BWA][bwa] by [Heng Li][lh3].

[bwa]: https://github.com/lh3/bwa
[gplv3]: http://www.gnu.org/licenses/gpl-3.0.html
[lh3]: https://github.com/lh3
