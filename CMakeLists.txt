cmake_minimum_required(VERSION 3.16)
project(montgomery_ops C ASM)

include(CheckCCompilerFlag)

if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message(FATAL_ERROR "montgomery_ops requires Linux/AArch64 with SVE2 support.")
endif()

if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)$")
  message(FATAL_ERROR "montgomery_ops requires an AArch64 target processor.")
endif()

check_c_compiler_flag("-march=armv8-a+sve2" MONTGOMERY_OPS_HAS_SVE2_FLAG)
if(NOT MONTGOMERY_OPS_HAS_SVE2_FLAG)
  message(FATAL_ERROR "The active compiler does not support -march=armv8-a+sve2.")
endif()

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

set(MONTGOMERY_OPS_ARCH_FLAGS "-march=armv8-a+sve2")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MONTGOMERY_OPS_ARCH_FLAGS} -Wall -Wextra -O2")
set(CMAKE_ASM_FLAGS "${CMAKE_ASM_FLAGS} ${MONTGOMERY_OPS_ARCH_FLAGS}")

include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src/include)

set(MONTGOMERY_OPS_SOURCES
    test.c
    src/poly_pointwise_montgomery_ref.c
    src/polyvecl_pointwise_acc_montgomery_l4_ref.c
    src/polyvecl_pointwise_acc_montgomery_l5_ref.c
    src/polyvecl_pointwise_acc_montgomery_l7_ref.c
    src/pointwise_montgomery.S
    src/mld_polyvecl_pointwise_acc_montgomery_l4.S
    src/mld_polyvecl_pointwise_acc_montgomery_l5.S
    src/mld_polyvecl_pointwise_acc_montgomery_l7.S
    src/pointwise_montgomery_sve2.S
    src/mld_polyvecl_pointwise_acc_montgomery_l4_sve2.S
    src/mld_polyvecl_pointwise_acc_montgomery_l5_sve2.S
    src/mld_polyvecl_pointwise_acc_montgomery_l7_sve2.S
)

add_executable(test_montgomery_ops ${MONTGOMERY_OPS_SOURCES})

add_custom_target(run
    COMMAND test_montgomery_ops
    DEPENDS test_montgomery_ops
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running test_montgomery_ops"
)
