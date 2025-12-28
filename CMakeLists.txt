cmake_minimum_required(VERSION 3.15)
project(stdcellgen LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(MSVC)
  add_compile_options(/W4 /permissive- /EHsc)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_library(stdcell
  src/core/core.cpp
  src/io/io.cpp
  src/matcher/matcher.cpp;src/matcher/matchdata.cpp;src/matcher/structural_matcher.cpp
  src/placer/placer.cpp
  src/router/router.cpp
)

target_include_directories(stdcell PUBLIC include)

add_executable(stdcellgen apps/stdcellgen.cpp)

target_link_libraries(stdcellgen PRIVATE stdcell)

# Install rules (optional)
install(TARGETS stdcellgen RUNTIME DESTINATION bin)
