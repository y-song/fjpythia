set(MODULE_NAME util)
set(MODULE ${PROJECT_NAME}_${MODULE_NAME})
set(NAME_LIB ${MODULE})

set(SOURCES_LIB
	strutil.cxx
	#sysutil.cxx
	looputil.cxx
	argparser.cxx
	fjutils.cxx
	pyutils.cxx
)

set_headers_lib(HEADERS_LIB SOURCES_LIB)

add_library(${NAME_LIB} SHARED ${SOURCES_LIB})
target_include_directories(${NAME_LIB} PUBLIC ${PYTHIA8_DIR}/include ${ROOT_INCLUDE_DIR} ${FASTJET_DIR}/include)
target_link_libraries(${NAME_LIB} ${PYTHIA8_LDFLAGS} ${ROOT_LIBRARIES} libRecursiveTools.a ${FASTJET_LIBS})

install_targets_custom()
