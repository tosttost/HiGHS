add_executable(highsexamplebuildmodel src/main.cpp)
target_link_libraries(highsexamplebuildmodel libhighs)

# install the binary
install(TARGETS highsexamplebuildmodel DESTINATION examples)