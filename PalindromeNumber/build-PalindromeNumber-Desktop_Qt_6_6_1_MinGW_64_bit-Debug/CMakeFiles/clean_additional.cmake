# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Debug")
  file(REMOVE_RECURSE
  "CMakeFiles\\PalindromeNumber_autogen.dir\\AutogenUsed.txt"
  "CMakeFiles\\PalindromeNumber_autogen.dir\\ParseCache.txt"
  "PalindromeNumber_autogen"
  )
endif()
