cmake_minimum_required(VERSION 3.12)
project(my-project)



list(APPEND CMAKE_PREFIX_PATH "D:/opencv/opencv/build")
find_package(OpenCV  REQUIRED)


include_directories(
	${OpenCV_INCLUDE_DIRS}
	${CMAKE_CURRENT_SOURCE_DIR}
)

include_directories($ENV{HOME}/.local/include/eigen3)

set(SOURCES
	PCA_FUNC.cpp
	PCA.h
	PCA_main.cpp
	main.cpp
)




add_executable(my-project ${SOURCES})


target_link_libraries(my-project ${OpenCV_LIBS})



