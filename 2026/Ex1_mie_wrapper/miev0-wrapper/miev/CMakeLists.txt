# SPDX-FileCopyrightText: 2026 Alberto P
#
# SPDX-License-Identifier: MPL-2.0

cmake_minimum_required(VERSION 3.16)
project(mie_scatt_bin LANGUAGES Fortran)

include(GNUInstallDirs)
enable_testing()

set(FORTRAN_SOURCES
    mieleg.F90
    MIEV0.f
    ErrPack.f
)
add_executable(miescat ${FORTRAN_SOURCES})

set_target_properties(miescat PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
)
target_compile_definitions(miescat PRIVATE
    $<$<CONFIG:Release>:NDEBUG>
)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(miescat PRIVATE
        $<$<CONFIG:Debug>:
            -O0                              # Disable optimizations for accurate debugging
            -g                               # Generate debugging information for GDB
            -fcheck=all                      # Enable all runtime checks (bounds, array, pointer, mem, do)
            -fbacktrace                      # Print a backtrace on runtime crash
            -ffpe-trap=invalid,zero,overflow # Trap on floating-point exceptions
            -finit-real=snan                 # Initialize REAL variables to Signaling NaN to catch uninitialized usage
            -finit-integer=-2147483648       # Initialize INTEGER variables to a large negative value (INT32_MIN)
            -Wall                            # Enable standard compiler warnings
            -Wextra                          # Enable extra compiler warnings
            -Wconversion                     # Warn about implicit type conversions
            -Wconversion-extra               # Warn about explicit/implicit conversions with precision or type changes
            -Wimplicit-interface             # Warn if a procedure is called without an explicit interface
        >
    )
endif()

install(TARGETS miescat
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

add_test(NAME test_mie
    COMMAND miescat 1.3484 0.001 1.0 0.41
)

find_program(VALGRIND_EXECUTABLE valgrind)
if(VALGRIND_EXECUTABLE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_test(NAME test_mie_valgrind
        COMMAND ${VALGRIND_EXECUTABLE}
            --leak-check=full
            --show-leak-kinds=all
            $<TARGET_FILE:miescat>
            1.3484 0.001 1.0 0.41
    )
endif()
