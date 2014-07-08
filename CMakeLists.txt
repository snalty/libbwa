cmake_minimum_required(VERSION 2.8)

project(bwa C)

set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
find_package(Threads)
find_package(ZLIB)

set(CMAKE_C_FLAGS "-g -Wall -O2 -Wno-unused-function")

set(BWA_HEADERS
  libbwa.h)

set(BWA_SOURCES
  libbwa.c
  libbwa_index.c
  libbwa_aln.c
  libbwa_samse.c
  libbwa_sampe.c
  libbwa_sw.c
  libbwa_mem.c
  utils.c
  kthread.c
  kstring.c
  ksw.c
  bwt.c
  bntseq.c
  bwa.c
  bwamem.c
  bwamem_pair.c
  bwamem_extra.c
  malloc_wrap.c
  QSufSort.c
  bwt_gen.c
  bwase.c
  bwaseqio.c
  bwtgap.c
  bwtaln.c
  bamlite.c
  is.c
  bwtindex.c
  bwape.c
  kopen.c
  pemerge.c
  bwtsw2_core.c
  bwtsw2_main.c
  bwtsw2_aux.c
  bwt_lite.c
  bwtsw2_chain.c
  fastmap.c
  bwtsw2_pair.c)

include_directories(${PROJECT_SOURCE_DIR})

set(SERIAL "0.7.9a")
set(SOSERIAL "0.7.9a")
add_library(bwa SHARED main.c ${BWA_SOURCES})
set_target_properties(bwa PROPERTIES VERSION ${SERIAL} SOVERSION ${SOSERIAL} MACOSX_RPATH 1)
target_link_libraries(bwa ${CMAKE_THREAD_LIBS_INIT} m z)

if(BUILD_TESTING)
  enable_testing()
  add_executable(bwa_test test.c ${BWA_SOURCES})
  add_test(bwa_test bwa_test)
  target_link_libraries(bwa_test ${CMAKE_THREAD_LIBS_INIT} m z cunit)
endif(BUILD_TESTING)

install(TARGETS bwa LIBRARY DESTINATION lib)
install(FILES ${BWA_HEADERS} DESTINATION include)
