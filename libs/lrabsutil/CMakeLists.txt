add_library(lrabsutil STATIC
    Assert.hh
    BoolVector.hh
    BoolVector.icc
    ConfigMap.hh
    ConfigMap.icc
    Exception.hh
    Exception.icc
    Hashes.hh
    Math.hh
    Math.icc
    Random.hh
    Random.icc
    Resources.hh
    SetIterators.hh
    SimpleStack.hh
    SimpleStack.icc
    String.hh
    String.icc
    UnionFind.hh
    ValidateableArray.hh
    ValidateableArray.icc
    Vector.hh
    Vector.icc
    BoolVector.cc
    ConfigMap.cc
    Exception.cc
    Random.cc
    Resources.cc
    String.cc
    UnionFind.cc)

get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
target_include_directories(lrabsutil PUBLIC ${PARENT_DIR})

install(TARGETS lrabsutil
        RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin
        LIBRARY DESTINATION ${CMAKE_BINARY_DIR}/bin
        ARCHIVE DESTINATION ${CMAKE_BINARY_DIR}/bin
       )

