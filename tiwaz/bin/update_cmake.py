#! /usr/bin/env python3

import glob
import os
import errno
import os.path


def find_files(folder, suffixes):
    files = sum((glob.glob("{}*{}".format(folder, s)) for s in suffixes), [])
    return list(sorted(files))


def find_source_files(folder):
    suffixes = [".c", ".C", ".cpp", ".c++"]
    return find_files(folder, suffixes)


def find_subdirectories(folder):
    dirs = sorted(glob.glob(folder + "*/"))
    return [d for d in dirs if "CMake" not in d]


def format_sources(target, sources):
    ret = ""

    line_pattern = "  PRIVATE ${{CMAKE_CURRENT_LIST_DIR}}/{:s}\n"

    if sources:
        ret += "".join(
            [
                "target_sources(" + target + "\n",
                "".join(line_pattern.format(os.path.basename(s)) for s in sources),
                ")\n\n",
            ]
        )

    return ret


def add_subdirectory(folder):
    line_pattern = "add_subdirectory({:s})\n"
    return line_pattern.format(os.path.basename(folder[:-1]))


def remove_file(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise


def append_to_file(filename, text):
    with open(filename, "a") as f:
        f.write(text)


def recurse(base_directory, targets):
    cmake_file = base_directory + "CMakeLists.txt"
    remove_file(cmake_file)

    source_files = find_source_files(base_directory)

    for dependency, target in targets.items():
        filtered_sources = list(filter(select_for(dependency), source_files))

        if dependency == "mpi":
            append_to_file(cmake_file, "if(ZISA_HAS_MPI)\n\n")

        append_to_file(cmake_file, format_sources(target, filtered_sources))

        if dependency == "mpi":
            append_to_file(cmake_file, "endif()\n")

    for d in find_subdirectories(base_directory):
        recurse(d, targets)
        append_to_file(cmake_file, add_subdirectory(base_directory + d))


def is_mpi_file(path):
    return "zisa/mpi/" in path


def is_generic_file(path):
    return not is_mpi_file(path)


def select_for(dependency):
    if dependency == "mpi":
        return is_mpi_file

    elif dependency == "generic":
        return is_generic_file

    else:
        raise Exception(f"Unknown dependency. [{dependency}]")


def add_executable(cmake_file, target, source_file):
    line = """
target_sources({}
  PRIVATE ${{CMAKE_CURRENT_SOURCE_DIR}}/{}
)
""".format(
        target, source_file
    )
    append_to_file(cmake_file, line)


if __name__ == "__main__":

    cmake_file = "src/CMakeLists.txt"
    remove_file(cmake_file)

    base_directory = "src/"
    for d in find_subdirectories(base_directory):
        recurse(d, {"generic": "zisa_generic_obj", "mpi": "zisa_mpi_obj"})
        append_to_file(cmake_file, add_subdirectory(base_directory + d))

    add_executable(cmake_file, "zisa", "zisa.cpp")
    add_executable(cmake_file, "locate-point", "locate_point.cpp")
    add_executable(cmake_file, "domain-decomposition", "domain_decomposition.cpp")
    add_executable(cmake_file, "opengl-demo", "opengl_demo.cpp")
    add_executable(cmake_file, "something-demo", "something_demo.cpp")
    add_executable(cmake_file, "renumber-grid", "renumber_grid.cpp")
    add_executable(
        cmake_file, "test-distributed-reference", "test_distributed_reference.cpp"
    )
    add_executable(
        cmake_file, "space-filling-curve-demo", "space_filling_curve_demo.cpp"
    )

    recurse("test/", {"generic": "unit_tests"})
    recurse("benchmarks/", {"generic": "micro_benchmarks"})
