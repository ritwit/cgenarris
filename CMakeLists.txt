cmake_minimum_required(VERSION  3.10.2)
project("cgenarris")

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_C_FLAGS "-Wall -Wpedantic -O3 -std=gnu99")

include_directories("${PROJECT_SOURCE_DIR}/src")
set(SOURCES
src/algebra.c
src/check_structure.c
src/combinatorics.c
src/crystal_utils.c
src/lattice_generator.c
src/molecule_placement.c
src/molecule_utils.c
src/pygenarris_mpi.c
src/randomgen.c
src/read_input.c
src/spg_generation.c

# Layer groups
src/lattice_generator_layer.c
src/find_layer_group.c

# spglib sources
src/spglib_src/arithmetic.c
src/spglib_src/cell.c
src/spglib_src/debug.c
src/spglib_src/delaunay.c
src/spglib_src/determination.c
src/spglib_src/hall_symbol.c
src/spglib_src/kgrid.c
src/spglib_src/kpoint.c
src/spglib_src/mathfunc.c
src/spglib_src/niggli.c
src/spglib_src/overlap.c
src/spglib_src/pointgroup.c
src/spglib_src/primitive.c
src/spglib_src/refinement.c
src/spglib_src/sitesym_database.c
src/spglib_src/site_symmetry.c
src/spglib_src/spacegroup.c
src/spglib_src/spg_database.c
src/spglib_src/spglib.c
src/spglib_src/spglib_f.c
src/spglib_src/spin.c
src/spglib_src/symmetry.c
src/spglib_src/tetrahedron_method.c
)

add_executable(${PROJECT_NAME} ${SOURCES} src/cgenarris_mpi.c)
add_executable(cgenarris_test ${SOURCES} src/test_main.c)
target_link_libraries(${PROJECT_NAME} m)
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})
target_link_libraries(cgenarris_test m)
target_link_libraries(cgenarris_test ${MPI_LIBRARIES})

enable_testing()
add_test(NAME cgenarris_test
         COMMAND cgenarris_test
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/regression_test1/)

