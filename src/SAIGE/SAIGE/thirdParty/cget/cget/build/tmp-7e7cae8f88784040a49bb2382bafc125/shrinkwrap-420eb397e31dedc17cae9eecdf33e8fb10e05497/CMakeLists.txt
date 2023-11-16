cmake_minimum_required(VERSION 3.2)
include(ExternalProject)
include(CMakePackageConfigHelpers)
project(shrinkwrap VERSION 1.0.0)

enable_testing()

set(CMAKE_CXX_STANDARD 11)

if (BUILD_SHARED_LIBS)
    set(LIBLZMA_LIB_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}lzma${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(ZLIB_LIB_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}z${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(ZSTD_LIB_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}zstd${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
    set(LIBLZMA_LIB_NAME ${CMAKE_STATIC_LIBRARY_PREFIX}lzma${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(ZLIB_LIB_NAME ${CMAKE_STATIC_LIBRARY_PREFIX}z${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(ZSTD_LIB_NAME ${CMAKE_STATIC_LIBRARY_PREFIX}zstd${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

find_library(LIBLZMA_LIBRARIES
             NAMES ${LIBLZMA_LIB_NAME})

find_library(ZLIB_LIBRARIES
             NAMES ${ZLIB_LIB_NAME})

find_library(ZSTD_LIBRARIES
             NAMES ${ZSTD_LIB_NAME})

if (NOT LIBLZMA_LIBRARIES)
    message(WARNING "lzma library not found")
endif()

if (NOT ZLIB_LIBRARIES)
    message(WARNING "zlib library not found")
endif()

if (NOT ZSTD_LIBRARIES)
    message(WARNING "zstd library not found")
endif()

add_library(shrinkwrap INTERFACE)
if (CMAKE_VERSION VERSION_GREATER 3.3)
    target_sources(shrinkwrap INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/shrinkwrap/stdio.hpp;${CMAKE_CURRENT_SOURCE_DIR}/include/shrinkwrap/xz.hpp;${CMAKE_CURRENT_SOURCE_DIR}/include/shrinkwrap/gz.hpp;${CMAKE_CURRENT_SOURCE_DIR}/include/shrinkwrap/zstd.hpp;${CMAKE_CURRENT_SOURCE_DIR}/include/shrinkwrap/istream.hpp>)
    target_include_directories(shrinkwrap INTERFACE
                               $<INSTALL_INTERFACE:include>
                               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
    target_link_libraries(shrinkwrap INTERFACE ${LIBLZMA_LIBRARIES} ${ZLIB_LIBRARIES} ${ZSTD_LIBRARIES})

    if (BUILD_TESTS)
        add_executable(shrinkwrap-test src/test.cpp)
        target_link_libraries(shrinkwrap-test shrinkwrap)
    endif()
else()
    if (BUILD_TESTS)
        add_executable(shrinkwrap-test src/test.cpp)
        target_link_libraries(shrinkwrap-test ${LIBLZMA_LIBRARIES} ${ZLIB_LIBRARIES} ${ZSTD_LIBRARIES})
        target_include_directories(shrinkwrap-test PUBLIC include)
    endif()
endif()

if (BUILD_TESTS)
    add_test(stdio_iterator_test shrinkwrap-test stdio-iter)
    add_test(stdio_seek_test shrinkwrap-test stdio-seek)
    add_test(xz_seek_test shrinkwrap-test xz-seek)
    add_test(xz_iterator_test shrinkwrap-test xz-iter)
    add_test(gz_iterator_test shrinkwrap-test gz-iter)
    add_test(bgzf_seek_test shrinkwrap-test bgzf-seek)
    add_test(bgzf_iterator_test shrinkwrap-test bgzf-iter)
    add_test(zstd_iterator_test shrinkwrap-test zstd-iter)
    add_test(zstd_seek_test shrinkwrap-test zstd-seek)
    add_test(generic_iterator_test shrinkwrap-test generic-iter)
    add_test(generic_seek_test shrinkwrap-test generic-seek)
endif()
install(DIRECTORY include/shrinkwrap DESTINATION include)
if (CMAKE_VERSION VERSION_GREATER 3.3)
    install(TARGETS shrinkwrap EXPORT ${PROJECT_NAME}-config
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)

    install(EXPORT ${PROJECT_NAME}-config DESTINATION share/${PROJECT_NAME})
    write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake COMPATIBILITY SameMajorVersion)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake DESTINATION share/${PROJECT_NAME})
    export(EXPORT ${PROJECT_NAME}-config)
endif()