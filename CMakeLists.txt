cmake_minimum_required(VERSION 3.5)
project(utilities)

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
if (WIN32)
	set(CMAKE_C_FLAGS "-D_WIN32")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fwrapv -ffast-math -Wall -Wextra -pedantic")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g3 -DDEBUG")

set(HEADERS
	Spawn.h
	Climates.h
	U_Math.h
)
set(SOURCES
	Spawn_Lookup_Tables.c
	Spawn.c
	Climates.c
	U_Math.c
)

add_subdirectory(cubiomes)

add_library(utilitiesObjects OBJECT ${SOURCES})
set_property(TARGET utilitiesObjects PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(utilities SHARED $<TARGET_OBJECTS:utilitiesObjects>)
add_library(utilities_static STATIC $<TARGET_OBJECTS:utilitiesObjects>)

target_link_libraries(utilities PUBLIC cubiomes)
target_include_directories(utilities PRIVATE HEADERS)

install(TARGETS utilities utilities_static DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)


