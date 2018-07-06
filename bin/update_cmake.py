#! /usr/bin/env python3

import glob
import os
import errno
import os.path

def find_files(folder, suffix):
    files = sorted(glob.glob(folder + "*" + suffix))
    return [os.path.basename(f) for f in files if "CMake" not in f]

def find_cpp_files(folder):
    return find_files(folder, ".cpp")

def find_cuda_files(folder):
    return find_files(folder, ".cu")

def find_subdirectories(folder):
    dirs = sorted(glob.glob(folder + "*/"))
    return [d for d in dirs if "CMake" not in d]

def format_sources(target, sources):
    ret = ""

    line_pattern = "  PRIVATE ${{CMAKE_CURRENT_LIST_DIR}}/{:s}\n"
    if sources:
        ret += "".join(["target_sources(" + target + "\n",
                        "".join(line_pattern.format(s) for s in sources),
                        ")\n\n"])

    return ret

def add_subdirectory(folder):
    line_pattern = "add_subdirectory({:s})\n"
    return  line_pattern.format(os.path.basename(folder[:-1]))

def truncate_file(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise

def append_to_file(filename, text):
    with open(filename, 'a') as f:
        f.write(text)

def recurse(base_directory, target):
    filename = base_directory + "CMakeLists.txt"
    write_output = lambda s : append_to_file(filename, s)
    truncate_file(filename)

    if target == "zisa":
        cpp = find_cpp_files(base_directory)
        write_output(format_sources(target, cpp))

        cuda = find_cuda_files(base_directory)
        if cuda:
            write_output("if(USE_CUDA)\n\n")
            write_output(format_sources(target, cuda))
            write_output("endif()\n\n")

    elif target == "unit_tests":
        cpp = find_cpp_files(base_directory)
        write_output(format_sources(target, cpp))

    else:
        raise Exception("Unknown 'target'. [{:s}]".format(target))

    subdirectories = find_subdirectories(base_directory)

    for subdirectory in subdirectories:
        write_output(add_subdirectory(base_directory + subdirectory))

    for subdirectory in subdirectories:
        recurse(subdirectory, target)

def add_executable(filename):
    line = """
target_sources(run
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/zisa.cpp
)
"""
    append_to_file(filename, line)


if __name__ == "__main__":

    filename = "src/CMakeLists.txt"
    base_directory = "src/"

    truncate_file(filename)
    subdirectories = find_subdirectories(base_directory)

    regular_folders = ["zisa/"]
    for d in regular_folders:
        recurse("src/" + d, "zisa")
        append_to_file(filename, add_subdirectory(base_directory + d))

    add_executable(filename)

    recurse("test/", "unit_tests")
