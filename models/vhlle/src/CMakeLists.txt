# Compile everything except the main source file into a static lib to be linked
# to both the main executable and the tests.
add_library(${LIBRARY_NAME} STATIC
fast.cpp  cll.cpp eos.cpp eos_hotqcd.cpp trancoeff.cpp fld.cpp hdo.cpp  main.cpp IC_reader.cpp rmn.cpp cornelius.cpp 
)

# Compile the actual executable.
set(MAIN main.cpp)
add_executable(${PROJECT_NAME} ${MAIN})
target_link_libraries(${PROJECT_NAME} ${LIBRARY_NAME} ${ROOT_LIBS} ${HDF5_LIBRARIES})

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
