include(../RibiLibraries/Apfloat.pri)
include(../RibiLibraries/BoostAll.pri)

include(plane.pri)
include(plane_apfloat.pri)
include(../RibiClasses/CppGeometry/CppGeometry.pri)
include(../RibiClasses/CppRibiRegex/CppRibiRegex.pri)

#
SOURCES += main_test_console.cpp

# C++14
QMAKE_CXX = g++-5
QMAKE_LINK = g++-5
QMAKE_CC = gcc-5
CONFIG += c++14
QMAKE_CXXFLAGS += -std=c++14

# High warning levels
# -Wshadow does not go with apfloat
QMAKE_CXXFLAGS += -Wall -Wextra -Wnon-virtual-dtor -pedantic -Werror


CONFIG(release, debug|release) {

  DEFINES += NDEBUG

  # gprof
  QMAKE_CXXFLAGS += -pg
  QMAKE_LFLAGS += -pg
}

CONFIG(debug, debug|release) {

  # gcov
  QMAKE_CXXFLAGS += -fprofile-arcs -ftest-coverage
  LIBS += -lgcov

  # UBSAN
  QMAKE_CXXFLAGS += -fsanitize=undefined
  QMAKE_LFLAGS += -fsanitize=undefined
  LIBS += -lubsan
}

# Boost.Test
LIBS += -lboost_unit_test_framework

# Boost.Graph
LIBS += \
  -lboost_date_time \
  -lboost_graph \
  -lboost_regex

# QResources give this error
QMAKE_CXXFLAGS += -Wno-unused-variable
