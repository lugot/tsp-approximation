SET(CPLEX_ROOT_DIR "" CACHE PATH "CPLEX root directory")

FIND_PATH(CPLEX_INCLUDE_DIR
    ilcplex/cplex .h
    PATHS "/opt/ibm/ILOG/CPLEX_Studio201/cplex/include"
    HINTS ${CPLEX_ROOT_DIR}/include
    )

FIND_LIBRARY(CPLEX_LIBRARY
    libcplex.a
    PATHS "/opt/ibm/ILOG/CPLEX_Studio201/cplex/lib/x86-64_linux/static_pic"
    HINTS ${CPLEX_ROOT_DIR}/lib)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CPLEX DEFAULT_MSG CPLEX_LIBRARY CPLEX_INCLUDE_DIR)

FIND_PATH(CPLEX_BIN_DIR
    cplex
    PATHS "/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux"
    )

IF (CPLEX_FOUND)
    SET(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR})
    SET(CPLEX_LIBRARIES ${CPLEX_LIBRARY})
    IF (CMAKE_SYSTEM_NAME STREQUAL "Linux")
        SET(CPLEX_LIBRARIES "${CPLEX_LIBRARIES};m;pthread")
    ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
ENDIF(CPLEX_FOUND)

MARK_AS_ADVANCED(CPLEX_LIBRARY CPLEX_INCLUDE_DIR CPLEX_BIN_DIR)
