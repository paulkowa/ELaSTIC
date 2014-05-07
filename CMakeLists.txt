CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(ELaSTIC)

SET(ELASTIC_VERSION_MAJOR 1)
SET(ELASTIC_VERSION_MINOR 65)

OPTION(WITH_SKETCH "Build elastic-sketch component" ON)
OPTION(WITH_TOOLS "Build all tools except of elastic-sketch" ON)

IF (CMAKE_CXX_FLAGS)
ELSE()
  SET(CMAKE_CXX_FLAGS "-std=c++0x -O3")
ENDIF()

FIND_PACKAGE(Boost 1.48.0)

IF (Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ELSE()
  MESSAGE(FATAL_ERROR "Boost header libraries not found!\nUse -DBOOST_ROOT=path to specify alternative location")
ENDIF()

IF (WITH_TOOLS)
  SET(Boost_USE_MULTITHREADED OFF)
  FIND_PACKAGE(Boost 1.48.0 COMPONENTS filesystem iostreams system)
  IF (Boost_FOUND)
  ELSE()
    MESSAGE(FATAL_ERROR "Boost binary libraries not found!\nMake sure Boost binary components are installed")
  ENDIF()
ENDIF()

IF (WITH_SKETCH)
  FIND_PACKAGE(MPI)
  IF (MPI_FOUND)
    MESSAGE(STATUS "Set MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")
  ELSE()
    MESSAGE(FATAL_ERROR "MPI library not found!\nMake sure MPI compiler is in your PATH")
  ENDIF()
ENDIF()

INCLUDE_DIRECTORIES(include src)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(doc)

IF (WITH_SKETCH AND WITH_TOOLS)
  MESSAGE(STATUS "All ELaSTIC tools enabled")
ELSEIF (WITH_SKETCH)
  MESSAGE(STATUS "Only elastic-sketch enabled")
ELSE()
  MESSAGE(STATUS "All ELaSTIC tools except of elastic-sketch enabled")
ENDIF()
