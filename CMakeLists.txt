cmake_minimum_required(VERSION 3.10)  # Minimum CMake version supporting C++17

project(scLJA
    VERSION 1.0
    LANGUAGES CXX
)

# Set C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # Don't use compiler-specific extensions

# Find OpenMP package
find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)

# Output directory for binaries (before installation)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Configure compiler warnings
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /WX")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Werror")
endif()

# Collect all source files from src directory
file(GLOB_RECURSE SOURCES "src/*.cpp")
file(GLOB_RECURSE HEADERS "src/*.hpp" "src/*.h")

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_link_libraries(scLJA PRIVATE OpenMP::OpenMP_CXX ${ZLIB_LIBRARIES})

# For older compilers that need -lstdc++fs
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(${PROJECT_NAME} PRIVATE stdc++fs)
endif()

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE src)
target_include_directories(${PROJECT_NAME} PRIVATE ${ZLIB_INCLUDE_DIRS})

# Add any required libraries here
# find_package(Boost REQUIRED COMPONENTS filesystem system)
# target_link_libraries(${PROJECT_NAME} PRIVATE Boost::filesystem Boost::system)

# Installation setup
install(TARGETS scLJA
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

# Create install target
install(DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
    DESTINATION bin
    FILES_MATCHING PATTERN "*"
)

# Add tests if you have them (optional)
# enable_testing()
# add_test(NAME YourTest COMMAND ${PROJECT_NAME} --test)