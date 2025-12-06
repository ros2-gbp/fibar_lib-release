function(add_clang_format_check target_name)
    find_program(CLANG_FORMAT_EXECUTABLE clang-format)
    if(NOT CLANG_FORMAT_EXECUTABLE)
        message(WARNING "clang-format not found. Skipping format check.")
        return()
    endif()

    file(GLOB_RECURSE CXX_SOURCE_FILES
        "${CMAKE_SOURCE_DIR}/test/*.h"
        "${CMAKE_SOURCE_DIR}/test/*.hpp"
        "${CMAKE_SOURCE_DIR}/test/*.c"
        "${CMAKE_SOURCE_DIR}/test/*.cpp"
        "${CMAKE_SOURCE_DIR}/include/*.h"
        "${CMAKE_SOURCE_DIR}/include/*.hpp"
    )

    add_custom_target(${target_name}
        COMMAND ${CMAKE_COMMAND} -E echo "Checking clang-format compliance..."
        COMMAND ${CLANG_FORMAT_EXECUTABLE} --style=file --dry-run ${CXX_SOURCE_FILES}
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    )
endfunction()