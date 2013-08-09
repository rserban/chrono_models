cmake_minimum_required(VERSION 2.8)
project(chrono_models)
include(SBELUtils.cmake)

SET(ENABLE_UNIT_GPU FALSE   CACHE   BOOL    "Turn on for GPU models")
SET(ENABLE_UNIT_MPI FALSE   CACHE   BOOL    "Turn on for GPU models")

# Additional packages
enable_openmp_support()

# User specified path to ChronoEngine
SET(CHRONO_INC       ""   CACHE PATH "Where is ChronoEngine Source located")

# Include directories based on user path to ChronoEngine
SET(INC_DIRS ${CHRONO_INC} "${CHRONO_INC}/collision" "${CHRONO_INC}/collision/bullet")

# User specified path to ChronoEngine Libraries
SET(CHRONO_LIB_PATH  ""   CACHE PATH "Where is CHRONOLIB located")

# Find required Libraries
FIND_LIBRARY(CHRONO_LIB NAMES ChronoEngine PATHS ${CHRONO_LIB_PATH} REQUIRED)

# Set unit-specific options
IF(ENABLE_UNIT_GPU)
    enable_cuda_support()
    LIST(APPEND INC_DIRS "${CUDA_TOOLKIT_ROOT_DIR}/include" "${CUDA_SDK_ROOT_DIR}/common/inc")

    FIND_LIBRARY(CHRONO_LIB_GPU NAMES ChronoEngine_GPU PATHS ${CHRONO_LIB_PATH} REQUIRED)
    FIND_LIBRARY(CHRONO_LIB_POSTPROCESS NAMES ChronoEngine_POSTPROCESS PATHS ${CHRONO_LIB_PATH} REQUIRED)
    FIND_LIBRARY(CHRONO_LIB_OPENGL NAMES ChronoEngine_OPENGL PATHS ${CHRONO_LIB_PATH} REQUIRED)

    SET(GPU_LIBS
        ${CHRONO_LIB_POSTPROCESS}
        ${CHRONO_LIB_GPU}
        ${CHRONO_LIB_OPENGL}
    )
ENDIF()

IF(ENABLE_UNIT_MPI)
    enable_mpi_support()

    FIND_LIBRARY(CHRONO_LIB_MPI NAMES ChronoEngine_MPI PATHS ${CHRONO_LIB_PATH} REQUIRED)

    SET(MPI_LIBS
        ${CHRONO_LIB_MPI}
    )
ENDIF()

include_directories(${INC_DIRS})

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ")

# Add executables
IF (ENABLE_UNIT_GPU)
    ADD_EXECUTABLE(ball_drop models/gpu/ball_drop.cpp)		
    TARGET_LINK_LIBRARIES(ball_drop	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(convergence models/gpu/convergence.cpp)		
    TARGET_LINK_LIBRARIES(convergence	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(aor models/gpu/aor.cpp)		
    TARGET_LINK_LIBRARIES(aor	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(diagonal_impact models/gpu/diagonal_impact.cpp)		
    TARGET_LINK_LIBRARIES(diagonal_impact	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(shaker models/gpu/shaker.cpp)		
    TARGET_LINK_LIBRARIES(shaker	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(character_walk models/gpu/character_walk.cpp)		
    TARGET_LINK_LIBRARIES(character_walk	${CHRONO_LIB} ${GPU_LIBS})

    ADD_EXECUTABLE(trough models/gpu/trough.cpp)		
    TARGET_LINK_LIBRARIES(trough	${CHRONO_LIB} ${GPU_LIBS})

    ADD_EXECUTABLE(density models/gpu/density.cpp)		
    TARGET_LINK_LIBRARIES(density	${CHRONO_LIB} ${GPU_LIBS})

    ADD_EXECUTABLE(fluid_tank models/gpu/fluid_tank.cpp)		
    TARGET_LINK_LIBRARIES(fluid_tank	${CHRONO_LIB} ${GPU_LIBS})

    ADD_EXECUTABLE(monkey models/gpu/monkey.cpp)		
    TARGET_LINK_LIBRARIES(monkey	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(sticky models/gpu/sticky.cpp)		
    TARGET_LINK_LIBRARIES(sticky	${CHRONO_LIB} ${GPU_LIBS})

    ADD_EXECUTABLE(foam models/gpu/foam.cpp)		
    TARGET_LINK_LIBRARIES(foam	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(rollers models/gpu/rollers.cpp)		
    TARGET_LINK_LIBRARIES(rollers	${CHRONO_LIB} ${GPU_LIBS})

    ADD_EXECUTABLE(sls models/gpu/sls.cpp)		
    TARGET_LINK_LIBRARIES(sls	${CHRONO_LIB} ${GPU_LIBS})
    
   # ADD_EXECUTABLE(walker models/generic/demo_walker_nonIrr.cpp)		
    #TARGET_LINK_LIBRARIES(walker	${CHRONO_LIB} ${GPU_LIBS})

    #ADD_EXECUTABLE(cylinder_drop models/gpu/cylinder_drop.cpp)		
    #TARGET_LINK_LIBRARIES(cylinder_drop	${CHRONO_LIB} ${GPU_LIBS})
ENDIF(ENABLE_UNIT_GPU)

IF (ENABLE_UNIT_MPI)
    ADD_EXECUTABLE(rover_enormous models/mpi/rover_enormous.cpp)
    TARGET_LINK_LIBRARIES(rover_enormous ${CHRONO_LIB} ${MPI_LIBS})
ENDIF()
