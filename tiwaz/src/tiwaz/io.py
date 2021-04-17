import json
import os
import pickle


def ensure_directory_exists(filename):
    dirname = os.path.dirname(filename)
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)


def read_something(filename, command, mode="r"):
    with open(filename, mode) as f:
        return command(f)


def write_something(filename, command, mode="w"):
    with open(filename, mode) as f:
        command(f)


def read_txt(filename):
    return read_something(filename, lambda f: f.read())


def write_txt(filename, string):
    write_something(filename, lambda f: f.write(string))


def read_json(filename):
    return read_something(filename, lambda f: json.load(f))


def write_json(filename, obj):
    write_something(filename, lambda f: json.dump(obj, f, indent=4))


def list_of_lists(ll):
    return [list(l) for l in ll]


def read_pickle(filename):
    return read_something(filename, lambda f: pickle.load(f), mode="rb")


def write_pickle(filename, obj):
    write_something(filename, lambda f: pickle.dump(obj, f), mode="wb")
