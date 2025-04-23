cmake_minimum_required(VERSION 3.10)
project(scLJA
    VERSION 1.0
    LANGUAGES CXX
)

#----------------------------------------------------------------
# 1. Basic Configuration
#----------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set output directories (all binaries go to bin/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Compiler warnings
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /WX")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Werror")
endif()

#----------------------------------------------------------------
# 2. Conda Environment Setup
#----------------------------------------------------------------
if(NOT DEFINED ENV{CONDA_PREFIX})
    message(FATAL_ERROR "Conda environment not activated! Please run:
    conda activate sclja")
else()
    message(STATUS "Using Conda environment: $ENV{CONDA_PREFIX}")
    
    # Set explicit Conda paths
    set(CONDA_LIB_DIR "$ENV{CONDA_PREFIX}/lib")
    set(CONDA_INCLUDE_DIR "$ENV{CONDA_PREFIX}/include")
    
    # Add to compiler search paths
    list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
    list(APPEND CMAKE_LIBRARY_PATH "${CONDA_LIB_DIR}")
    list(APPEND CMAKE_INCLUDE_PATH "${CONDA_INCLUDE_DIR}")
    
    # Ensure linker uses Conda libraries first
    link_directories("${CONDA_LIB_DIR}")
endif()

#----------------------------------------------------------------
# 3. Dependency Detection (Conda-only)
#----------------------------------------------------------------
find_package(OpenMP REQUIRED)

# Find all libraries strictly in Conda environment
find_library(HTS_LIB hts PATHS "${CONDA_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)
find_library(ZLIB_LIB z PATHS "${CONDA_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)
find_library(BZ2_LIB bz2 PATHS "${CONDA_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)
find_library(LZMA_LIB lzma PATHS "${CONDA_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)
find_library(CURL_LIB curl PATHS "${CONDA_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)
find_library(CRYPTO_LIB crypto PATHS "${CONDA_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)

# Find headers
find_path(HTS_INCLUDE_DIR hts.h PATHS "${CONDA_INCLUDE_DIR}/htslib" NO_DEFAULT_PATH REQUIRED)
find_path(ZLIB_INCLUDE_DIR zlib.h PATHS "${CONDA_INCLUDE_DIR}" NO_DEFAULT_PATH REQUIRED)

#----------------------------------------------------------------
# 4. Target Configuration
#----------------------------------------------------------------
file(GLOB_RECURSE SOURCES "src/*.cpp")
file(GLOB_RECURSE HEADERS "src/*.hpp" "src/*.h")

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

target_include_directories(${PROJECT_NAME} PRIVATE
    src
    ${HTS_INCLUDE_DIR}
    ${ZLIB_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    OpenMP::OpenMP_CXX
    ${HTS_LIB}
    ${ZLIB_LIB}
    ${BZ2_LIB}
    ${LZMA_LIB}
    ${CURL_LIB}
    ${CRYPTO_LIB}
)

# Handle older GCC versions
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(${PROJECT_NAME} PRIVATE stdc++fs)
endif()

#----------------------------------------------------------------
# 5. Installation (optional)
#----------------------------------------------------------------
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# Copy all binaries to project root bin/ (optional)
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_NAME}
    ${CMAKE_SOURCE_DIR}/bin/${PROJECT_NAME}
)