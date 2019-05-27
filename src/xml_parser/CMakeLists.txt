# Project setup
project(xml C CXX)
set(VERSION_MAJOR "0")
set(VERSION_MINOR "2")
set(VERSION_PATCH "0")
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR) 


# Define main library target
add_library(yxml STATIC "")

# Sources
target_sources(
	yxml
	PRIVATE
		"yxml.c"
)


target_include_directories(
	yxml
	PUBLIC
		"."
)

set_property(TARGET yxml PROPERTY POSITION_INDEPENDENT_CODE ON)
