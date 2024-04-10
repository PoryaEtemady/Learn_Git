# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Debug")
  file(REMOVE_RECURSE
  "CMakeFiles\\MergeSortedArray_autogen.dir\\AutogenUsed.txt"
  "CMakeFiles\\MergeSortedArray_autogen.dir\\ParseCache.txt"
  "MergeSortedArray_autogen"
  )
endif()
