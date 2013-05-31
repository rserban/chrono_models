cmake_minimum_required(VERSION 2.8)
project(hammad_models)
include(SBELUtils.cmake)

# Additional packages
enable_cuda_support()
enable_openmp_support()

# User specified path to ChronoEngine
SET(CHRONO_INC       ""   CACHE PATH "Where is ChronoEngine Source located")
# Include directories based on user path to ChronoEngine
include_directories( ${CHRONO_INC} "${CHRONO_INC}/collision" "${CHRONO_INC}/collision/bullet" "${CUDA_TOOLKIT_ROOT_DIR}/include" "${CUDA_SDK_ROOT_DIR}/common/inc" )
# User specified path to ChronoEngine Libraries
SET(CHRONO_LIB_PATH  ""   CACHE PATH "Where is CHRONOlIB located")
# Find required Libraries
FIND_LIBRARY(CHRONO_LIB NAMES ChronoEngine PATHS ${CHRONO_LIB_PATH} REQUIRED)
FIND_LIBRARY(CHRONO_LIB_GPU NAMES ChronoEngine_GPU PATHS ${CHRONO_LIB_PATH} REQUIRED)
FIND_LIBRARY(CHRONO_LIB_POSTPROCESS NAMES ChronoEngine_POSTPROCESS PATHS ${CHRONO_LIB_PATH} REQUIRED)
FIND_LIBRARY(CHRONO_LIB_OPENGL NAMES ChronoEngine_OPENGL PATHS ${CHRONO_LIB_PATH} REQUIRED)

set (LINK_LIBS
	${CHRONO_LIB}
	${CHRONO_LIB_POSTPROCESS}
	${CHRONO_LIB_GPU}
	${CHRONO_LIB_OPENGL}
)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ")

# Add executables
ADD_EXECUTABLE(ball_drop models/gpu/ball_drop.cpp)		
TARGET_LINK_LIBRARIES(ball_drop	${LINK_LIBS})

ADD_EXECUTABLE(convergence models/gpu/convergence.cpp)		
TARGET_LINK_LIBRARIES(convergence	${LINK_LIBS})

ADD_EXECUTABLE(aor models/gpu/aor.cpp)		
TARGET_LINK_LIBRARIES(aor	${LINK_LIBS})

ADD_EXECUTABLE(diagonal_impact models/gpu/diagonal_impact.cpp)		
TARGET_LINK_LIBRARIES(diagonal_impact	${LINK_LIBS})

ADD_EXECUTABLE(shaker models/gpu/shaker.cpp)		
TARGET_LINK_LIBRARIES(shaker	${LINK_LIBS})

#ADD_EXECUTABLE(cylinder_drop models/gpu/cylinder_drop.cpp)		
#TARGET_LINK_LIBRARIES(cylinder_drop	${LINK_LIBS})