
#---Adding example B2 subdirectories explicitly 
#   and a custom target to for building all example B2 options ----------

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

add_subdirectory(B3a)
add_subdirectory(B3b)
add_custom_target(B3 DEPENDS exampleB3a exampleB3b)

set( target_name B3 )
set(CMAKE_CXX_FLAGS "-g -O0")

message([STATUS] "Now checking for TACC Build ...")

if(DEFINED ENV{TACC})
	# TACC specific
	message([STATUS] "ON TACC")
	target_compile_definitions(B3 TACC=TRUE)
	# set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS $<TACC>)
endif()