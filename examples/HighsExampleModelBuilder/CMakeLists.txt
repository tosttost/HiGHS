add_executable(highsexamplemodelbuilder src/main.cpp)
target_link_libraries(highsexamplemodelbuilder libhighs)

# install the binary
install(TARGETS highsexamplemodelbuilder DESTINATION examples)