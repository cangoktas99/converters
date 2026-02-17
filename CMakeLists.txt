cmake_minimum_required(VERSION 3.10)
project(OF_GMF_Converters)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(LIBMESHB_DIR ${CMAKE_SOURCE_DIR}/../libMeshb-7.80/sources)
set(LIBMESHB_SRC ${LIBMESHB_DIR}/libmeshb7.c)

include_directories(${LIBMESHB_DIR})

add_executable(of2gmf of2gmf.cpp ${LIBMESHB_SRC})
target_link_libraries(of2gmf m)

add_executable(gmf2of gmf2of.cpp ${LIBMESHB_SRC})
target_link_libraries(gmf2of m)
