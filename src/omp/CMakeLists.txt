SET(PREPARE_PROG elastic-prepare)
SET(SKETCH_OMP_PROG elastic-sketch-omp)

FIND_PACKAGE(OpenMP)

IF (OPENMP_FOUND)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

ADD_EXECUTABLE(${PREPARE_PROG} elastic-prepare.cpp)
TARGET_LINK_LIBRARIES(${PREPARE_PROG} ${IOMANIP_LIB} ${Boost_LIBRARIES})

INCLUDE_DIRECTORIES(../mpi)

ADD_EXECUTABLE(${SKETCH_OMP_PROG} elastic-sketch-omp.cpp)
TARGET_LINK_LIBRARIES(${SKETCH_OMP_PROG} ${IOMANIP_LIB})

INSTALL(TARGETS ${PREPARE_PROG} RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
INSTALL(TARGETS ${PREPARE_OMP_PROG} RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
