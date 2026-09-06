# SPDX-FileCopyrightText: 2026 Alberto P
#
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.16)
project(mie_code LANGUAGES CXX Fortran)

find_package(Python REQUIRED COMPONENTS Interpreter)
include(GNUInstallDirs)
enable_testing()

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

add_subdirectory(miev)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(WIN32)
        target_link_options(miescat PRIVATE
            -static
            #-no-pthread
            #"-l:libwinpthread.a"
            #-static-libgfortran
            #-static-libgcc
            #-static-libquadmath
        )
    elseif(APPLE)
        target_link_options(miescat PRIVATE
            -static-libgfortran
            -static-libgcc
            -static-libquadmath
        )
    else()
        target_link_options(miescat PRIVATE
            -static-libgfortran
        )
    endif()
endif()

add_custom_target(python_package
    COMMAND ${Python_EXECUTABLE} -m pip wheel ${CMAKE_SOURCE_DIR} -w ${CMAKE_BINARY_DIR}/dist
    DEPENDS miescat
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Build Python wheel package via setuptools"
)

install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/scripts/mieleg.py
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    RENAME mieleg
)

add_test(NAME test_mie_py
    COMMAND Python::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/scripts/mieleg.py
       --executable ${CMAKE_CURRENT_BINARY_DIR}/miev/miescat
       --m_real 1.3484 --m_img 0.001 --radius 1.0 --wavelength 0.41
)
