cmake_minimum_required(VERSION 3.2)
project(xz VERSION 5.2.3)

execute_process(COMMAND ./configure --prefix=${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(xz ALL COMMAND make WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMENT "Building xz ...")

install(DIRECTORY src/liblzma/api/ DESTINATION include FILES_MATCHING PATTERN "*.h")

if (BUILD_SHARED_LIBS)
    install(FILES
            src/liblzma/.libs/${CMAKE_SHARED_LIBRARY_PREFIX}lzma${CMAKE_SHARED_LIBRARY_SUFFIX}
            src/liblzma/.libs/${CMAKE_SHARED_LIBRARY_PREFIX}lzma.5${CMAKE_SHARED_LIBRARY_SUFFIX}
            DESTINATION lib OPTIONAL)
else()
    install(FILES
            src/liblzma/.libs/${CMAKE_STATIC_LIBRARY_PREFIX}lzma${CMAKE_STATIC_LIBRARY_SUFFIX}
            DESTINATION lib)
endif()