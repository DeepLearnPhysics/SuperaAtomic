cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# general configuration for compilation:
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -O3 -g")


project(supera)

set(CMAKE_PACKAGE_DIR "./")
set(PROJECT_SOURCE_DIR "./")
set(PYBIND11_MKDOC_DOCSTRING_FILE docstrings.h)

include_directories("./")

message("    Configuring supera build")


#########
# Python
#########
if(WITH_PYBIND OR WITH_PYROOT)

    find_package(PythonLibs REQUIRED)
    message("    Adding Python")
    message("    INCLUDE: " , ${PYTHON_INCLUDE_DIRS})
    include_directories(${PYTHON_INCLUDE_DIRS})

    if(NOT PYTHON_EXECUTABLE)
        message( FATAL_ERROR "PYTHON_EXECUTABLE not found. Cannot build edep2supera...")
    endif()

    ########
    # Numpy
    ########
    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}" -c
        "from __future__ import print_function\ntry: import numpy; print(numpy.get_include(), end='')\nexcept:pass\n"
        OUTPUT_VARIABLE __numpy_path)

    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}" -c
        "from __future__ import print_function\ntry: import numpy; print(numpy.__version__, end='')\nexcept:pass\n"
        OUTPUT_VARIABLE __numpy_version)

    find_path(PYTHON_NUMPY_INCLUDE_DIR numpy/arrayobject.h
        HINTS "${__numpy_path}" "${PYTHON_INCLUDE_PATH}" NO_DEFAULT_PATH)

    if(PYTHON_NUMPY_INCLUDE_DIR)
        set(PYTHON_NUMPY_FOUND 1 CACHE INTERNAL "Python numpy found")
    else()
        message( FATAL_ERROR "Numpy could not be located. Cannot build edep2supera...")
    endif()

    message("    Adding Numpy")
    message("    INCLUDE: " ${PYTHON_NUMPY_INCLUDE_DIR})
    include_directories(${PYTHON_NUMPY_INCLUDE_DIR})

    if(WITH_PYROOT)

        ################
        # Add PyROOT
        ################
        find_package(ROOT REQUIRED
          COMPONENTS Geom Physics Matrix MathCore Tree RIO)
        if(ROOT_FOUND)
          include(${ROOT_USE_FILE})
        else()
          message( FATAL_ERROR "ROOT is needed but not found. Cannot build edep2supera...")
        endif(ROOT_FOUND)

    endif()

    if(WITH_PYBIND)
        ####################################
        # Add pybind11
        ####################################
        add_subdirectory(pybind11)
        message("    Adding pybind11")
        message("    INCLUDE: " ${PYBIND11_INCLUDE_DIR})
        include_directories(${PYBIND11_INCLUDE_DIR})

        ####################################
        # Check if pybind11_mkdoc is around
        ####################################
        execute_process(
                COMMAND "${Python_EXECUTABLE}" -c "import pybind11_mkdoc"
                RESULT_VARIABLE PYBIND11_MKDOC_EXIT_CODE
        )
        if (${PYBIND11_MKDOC_EXIT_CODE} EQUAL 0)
            message("Found pybind11_mkdoc.  Will use it to autogenerate Python docstrings.")
            set(USE_PYBIND11_MKDOC ON)
            add_definitions(-DUSE_MKDOC_DOCSTRINGS)
            include_directories("${CMAKE_INSTALL_PREFIX}")
        endif()

    endif()

endif()

#
# Go build supera
#
add_subdirectory(supera)
