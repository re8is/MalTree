# 1 -1 3 4 -1 -1  7 -1 -1
cmake_minimum_required(VERSION 3.0.0)

project(TEVL)

set(CMAKE_CXX_STANDARD 17)

set(HOME .)

set(SRC_LIST online_computation_cost.cpp Cipher.cpp Utils.cpp PRGTest.cpp)

# set(OPENSSL_ROOT_DIR "/usr/bin/openssl")

find_package(OpenSSL REQUIRED)

set(EXECUTABLE_OUTPUT_PATH ${HOME}/bin)




add_executable(app ${SRC_LIST})
message(${PROJECT_SOURCE_DIR})
target_include_directories(app PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(app PRIVATE OpenSSL::SSL OpenSSL::Crypto)