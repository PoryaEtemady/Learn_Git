# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Debug")
  file(REMOVE_RECURSE
  "CMakeFiles\\RemoveDuplicatesFromSortedArray2_autogen.dir\\AutogenUsed.txt"
  "CMakeFiles\\RemoveDuplicatesFromSortedArray2_autogen.dir\\ParseCache.txt"
  "RemoveDuplicatesFromSortedArray2_autogen"
  )
endif()
