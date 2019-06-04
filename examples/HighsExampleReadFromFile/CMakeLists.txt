add_executable(highsexamplereadfromfile src/main.cpp)
target_link_libraries(highsexamplereadfromfile libhighs)

# install the binary
install(TARGETS highsexamplereadfromfile DESTINATION examples)