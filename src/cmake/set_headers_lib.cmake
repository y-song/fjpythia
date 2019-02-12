FUNCTION(PREPEND var prefix)
   SET(listVar "")
   FOREACH(f ${ARGN})
      LIST(APPEND listVar "${prefix}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

FUNCTION(SET_HEADERS_LIB var)
   # Headers from sources
   #string(REPLACE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/../include HEADERS_LIB "${HEADERS_LIB1}")
   string(REPLACE ".cxx" ".h" L_HEADERS_LIB1 "${SOURCES_LIB}")
   # message(STATUS "CURRENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}")
   # message(STATUS "CMAKE SOURCE_DIR ${CMAKE_SOURCE_DIR}")
   string(REPLACE "${CMAKE_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}/../include/${PROJECT_NAME}" L_HEADERS_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
   PREPEND(L_HEADERS_LIB ${L_HEADERS_DIR} ${L_HEADERS_LIB1})
   # message(STATUS "L_HEADERS_DIR=${L_HEADERS_DIR}")
   # message(STATUS "L_HEADERS_LIB=${L_HEADERS_LIB}")
   SET(${var} "${L_HEADERS_LIB}" PARENT_SCOPE)
ENDFUNCTION(SET_HEADERS_LIB)