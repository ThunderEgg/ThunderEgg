# External Project Modules
This directory includes all the modules need to use ThunderEgg with CMake's 
[external project](https://cmake.org/cmake/help/latest/module/ExternalProject.html) 
feature.

`ThunderEggExternalProject.cmake` contains a functions for easily adding ThunderEgg as an external project.

Example use case:

	include(ThunderEggExternalProject)
	ThunderEggExternalProject(TAG v1.0.1 COMPONENTS P4EST)

This makes available the target "ThunderEgg::ThunderEgg" with p4est optional features enabled. See `ThunderEggExternalProject.cmake` for more documentation.
