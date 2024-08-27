set(name data)

#
# List source/headers
#
# Get all the source/header files:
file(GLOB SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cxx)
file(GLOB HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h)
list(REMOVE_ITEM HEADERS "LinkDef.h")

# update the list inherited from src/supera
list(APPEND ALL_HEADERS ${HEADERS})
set(ALL_HEADERS "${ALL_HEADERS}" PARENT_SCOPE)

#
# Python binding preparation
#
if(NOT WITH_PYBIND)

    list(FILTER SOURCES EXCLUDE REGEX ".*_pybind\\.cxx$")
    list(FILTER HEADERS EXCLUDE REGEX ".*_pybind\\.h$")

else()

    if(USE_PYBIND11_MKDOC)
        add_dependencies(${name} mkdoc_docstring_file)
    endif()

endif()

if(WITH_PYROOT)

    # Generate dictionary
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    ROOT_GENERATE_DICTIONARY(G__supera_${name}
        ${HEADERS}
        OPTIONS -inlineInputHeader
        LINKDEF LinkDef.h)
    list(APPEND SOURCES G__supera_${name}.cxx)

endif()

#
# Compile
#
add_library(${name} OBJECT ${SOURCES})


#
# Install
#

# Header files
install (FILES ${HEADERS}
    DESTINATION ${CMAKE_PACKAGE_DIR}/include/supera/${name})


if(WITH_PYROOT)
    # Rootmap (for PyROOT dictionary)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libsupera_${name}.rootmap DESTINATION ${CMAKE_PACKAGE_DIR}/lib/supera)

    # pcm (for PyROOT dictionary)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libsupera_${name}_rdict.pcm DESTINATION ${CMAKE_PACKAGE_DIR}/lib/supera)
endif()
