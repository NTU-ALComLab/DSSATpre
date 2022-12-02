cmake_minimum_required(VERSION 3.11)

project(DSSATpre VERSION 1.0 LANGUAGES C CXX)

include(CheckIncludeFiles)
include(CheckIPOSupported)
include(CheckFunctionExists)
include(CheckLibraryExists)
include(CheckTypeSize)
include(TestBigEndian)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")


# Enforce C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANARD_REQUIRED ON)
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()

# Enable warnings and debug information
set(CMAKE_EXE_LINKER_FLAGS "-Wall -Wextra -Wno-sign-conversion -g ${CMAKE_EXE_LINKER_FLAGS}" )
set(CMAKE_C_FLAGS          "-Wall -Wextra -Wno-sign-conversion -g ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS        "-Wall -Wextra -Wno-sign-conversion -g ${CMAKE_CXX_FLAGS}")


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Definitions for using Easylogging++
option(NO_LOGGING:BOOL "Disable logging" OFF)
if (NO_LOGGING)
    add_definitions("-DELPP_DISABLE_LOGS -Wno-unused-parameter")
    message("-- Logging disabled")
endif()

option(STATIC:BOOL "Use static linking" OFF)
if (STATIC)
    if (WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES}")
    else()
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    endif()
    message("-- Enabled static linking.")
endif()


# Check if link-time optimizations are supported
check_ipo_supported(RESULT USE_IPO OUTPUT error)


# Find libz (used by boost-iostreams)
find_package(ZLIB REQUIRED)
add_library(zlib IMPORTED INTERFACE)
set_property(TARGET zlib PROPERTY INTERFACE_LINK_LIBRARIES  "${ZLIB_LIBRARIES}")
set_property(TARGET zlib PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${ZLIB_INCLUDES}")


# Find librt (used for timeout handling in hqspre)
find_package(RT REQUIRED)
add_library(rt IMPORTED INTERFACE)
set_property(TARGET rt PROPERTY INTERFACE_LINK_LIBRARIES  "${RT_LIBRARIES}")
set_property(TARGET rt PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${RT_INCLUDES}")

find_package(LTDL REQUIRED)
add_library(ltdl IMPORTED INTERFACE)
set_property(TARGET ltdl PROPERTY INTERFACE_LINK_LIBRARIES  "${LTDL_LIBRARIES}")
set_property(TARGET ltdl PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${LTDL_INCLUDES}")


# Find Boost libraries
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
if (STATIC)
    set(Boost_USE_STATIC_LIBS ON)
else()
    set(Boost_USE_STATIC_LIBS OFF)
endif()
find_package(Boost 1.40 REQUIRED COMPONENTS iostreams program_options)
add_library(boost IMPORTED INTERFACE)
set_property(TARGET boost PROPERTY INTERFACE_LINK_LIBRARIES  "${Boost_LIBRARIES}" zlib)
set_property(TARGET boost PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}")


# Try to find GLPK (used as MILP solver)
find_package(GLPK REQUIRED)
add_library(glpk IMPORTED INTERFACE)
set_property(TARGET glpk PROPERTY INTERFACE_LINK_LIBRARIES      "${GLPK_LIBRARIES};${LTDL_LIBRARIES}")
set_property(TARGET glpk PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${GLPK_INCLUDE_DIRS}")


# Enable OpenMP (used by the Antom SAT solver)
find_package(OpenMP REQUIRED)
add_library(openmp IMPORTED INTERFACE)
set_property(TARGET openmp PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
set_property(TARGET openmp PROPERTY INTERFACE_LINK_LIBRARIES  ${OpenMP_CXX_FLAGS})

# Check if valgrind is installed. This is used by the macro
# val_assert(...) to print a stack trace when the condition
# is violated and run through valgrind.
check_include_files("valgrind/valgrind.h" HAVE_VALGRIND)
if(HAVE_VALGRIND)
  add_definitions("-DHAVE_VALGRIND")
endif()




# Check if Doxygen is available; if yes, add the build target 'doc'
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/src/Doxyfile @ONLY)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/src/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

file(GLOB CHECK_CXX_SOURCE_FILES ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/src/*.hpp
    ${PROJECT_SOURCE_DIR}/src/*.ipp
    ${PROJECT_SOURCE_DIR}/libs/aigpp/*.cc
    ${PROJECT_SOURCE_DIR}/libs/aigpp/*.hh
    ${PROJECT_SOURCE_DIR}/libs/aigpp/*.icc
)

find_program(CLANG_FORMAT "clang-format")
if(CLANG_FORMAT)
    add_custom_target(
        format
        COMMAND
        ${CLANG_FORMAT}
        -i
        -style=file
        ${CHECK_CXX_SOURCE_FILES}
        COMMENT "Auto formatting of all source files"
    )

    add_custom_target(
        check-format
        COMMAND
        ${CLANG_FORMAT}
        -style=file
        -output-replacements-xml
        ${CHECK_CXX_SOURCE_FILES}
        # print output
        | tee ${CMAKE_BINARY_DIR}/check_format_file.txt | grep -c "replacement " |
                tr -d "[:cntrl:]" && echo " replacements necessary"
        # WARNING: fix to stop with error if there are problems
        COMMAND ! grep -c "replacement "
                  ${CMAKE_BINARY_DIR}/check_format_file.txt > /dev/null
        COMMENT "Checking format compliance"
    )
endif()

check_include_files("assert.h" HAVE_ASSERT_H)
check_include_files("float.h" HAVE_FLOAT_H)
check_include_files("dlfcn.h" HAVE_DLFCN_H)
check_include_files("inttypes.h" HAVE_INTTYPES_H)
check_include_files("limits.h" HAVE_LIMITS_H)
check_include_files("math.h" HAVE_MATH_H)
check_include_files("memory.h" HAVE_MEMORY_H)
check_include_files("stddef.h" HAVE_STDDEF_H)
check_include_files("stdint.h" HAVE_STDINT_H)
check_include_files("stdlib.h" HAVE_STDLIB_H)
check_include_files("strings.h" HAVE_STRINGS_H)
check_include_files("string.h" HAVE_STRING_H)
check_include_files("unistd.h" HAVE_UNISTD_H)
check_include_files("sys/resource.h" HAVE_SYS_RESOURCE_H)
check_include_files("sys/stat.h" HAVE_SYS_STAT_H)
check_include_files("sys/times.h" HAVE_SYS_TIMES_H)
check_include_files("sys/time.h" HAVE_SYS_TIME_H)
check_include_files("sys/types.h" HAVE_SYS_TYPES_H)
check_include_files("sys/wait.h" HAVE_SYS_WAIT_H)
check_function_exists("gethostname" HAVE_GETHOSTNAME)
check_function_exists("getrlimit" HAVE_GETRLIMIT)
check_function_exists("getrusage" HAVE_GETRUSAGE)
check_library_exists("m" "pow" "math.h" HAVE_POW)
check_library_exists("m" "powl" "math.h" HAVE_POWL)
check_library_exists("m" "sqrt" "math.h" HAVE_SQRT)
check_function_exists("strchr" HAVE_STRCHR)
check_function_exists("strstr" HAVE_STRSTR)
check_function_exists("sysconf" HAVE_SYSCONF)
check_type_size("int" SIZEOF_INT)
check_type_size("long" SIZEOF_LONG)
check_type_size("long double" SIZEOF_LONG_DOUBLE)
check_type_size("void*" SIZEOF_VOID_P)
set(CMAKE_EXTRA_INCLUDE_FILES "stddef.h")
check_type_size("ptrdiff_t" PTRDIFF_T)
set(CMAKE_EXTRA_INCLUDE_FILES)
test_big_endian(WORDS_BIGENDIAN)
set(CUDD_VERSION "3.0.0")

configure_file(libs/cudd-3.0.0/cudd_config.h.in libs/cudd-3.0.0/cudd_config.h @ONLY)

add_subdirectory("libs")
add_subdirectory("src")
