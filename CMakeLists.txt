cmake_minimum_required(VERSION 3.12)
project(IntegratedKMC)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wno-unused-variable -Wno-unused-but-set-variable")

find_package(Boost REQUIRED COMPONENTS system)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${Boost_INCLUDE_DIRS})

# 关键：添加新的源文件
set(SOURCES
        src/main.cpp
        src/globals.cpp
        src/efficient_explicit_sequence_record_number.cpp
        src/depolymerization_explicit_sequence_record_number.cpp
        src/LengthDependentTerminationRate.cpp
        src/my_rand.cpp
        )

set(HEADERS
        include/DataStructureAndConstant.h
        include/efficient_explicit_sequence_record_number.h
        include/depolymerization_explicit_sequence_record_number.h
        include/GlassEffectPropagationRate.h
        include/InitiatorEfficiency.h
        include/LengthDependentTerminationRate.h
        include/MolecularWeightDistribution.h
        include/MolecularWeightOutput.h
        include/my_rand.h
        include/ode.h
        include/PolymerFeathure.h
        include/ReadPolymer.h
        )

add_executable(integrated_kmc ${SOURCES} ${HEADERS})
target_link_libraries(integrated_kmc ${Boost_LIBRARIES})

target_compile_definitions(integrated_kmc PRIVATE
        BOOST_SYSTEM_NO_DEPRECATED
        _USE_MATH_DEFINES
        BOOST_ALLOW_DEPRECATED_HEADERS
        )

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output/poly_)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output/depoly_)

configure_file(${CMAKE_SOURCE_DIR}/src/input0.txt ${CMAKE_BINARY_DIR}/src/input0.txt COPYONLY)