cmake_minimum_required(VERSION 3.9)
project(Coursach)

set(CMAKE_CXX_STANDARD 11)
set(LIBRARIES
        Calc/lib/euler.cpp
        Calc/lib/fmin.cpp
        Calc/lib/lagrange.cpp
        Calc/lib/matrix.cpp
        Calc/lib/quanc8.cpp
        Calc/lib/rk3.cpp
        Calc/lib/rkf45.cpp
        Calc/lib/spline.cpp
        Calc/lib/ZEROIN.cpp)


add_library(lib STATIC ${LIBRARIES})


add_executable(firstLab ../FirstLab/FirstLab.cpp)
add_executable(secondLab ../SecondLab/SecondLab.cpp)
add_executable(thirdLab ../ThirdLab/ThirdLab.cpp)
add_executable(coursach ../Coursach/Coursach/Coursach.cpp)


target_link_libraries(firstLab lib)
target_link_libraries(secondLab lib)
target_link_libraries(thirdLab lib)
target_link_libraries(coursach lib)
