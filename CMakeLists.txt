project(my-project)


find_package(OpenCV REQUIRED)


include_directories(
	${OpenCV_INCLUDE_DIRS}
	${CMAKE_CURRENT_SOURCE_DIR}
)



set(SOURCES
	PCA_FUNC.cpp
	PCA.h
	PCA_main.cpp
	main.cpp
)




add_executable(my-project ${SOURCES})


target_link_libraries(my-project ${OpenCV_LIBS})



