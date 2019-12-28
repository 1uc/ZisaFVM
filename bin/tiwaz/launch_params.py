import itertools
import os
import subprocess
import glob
import multiprocessing

from .site_details import get_host
from .utils import merge_dict


class LaunchParams(object):
    """Iterable of launch parameters of zisa."""

    def __init__(self, params=None):
        self.params = params

    def __iter__(self):
        return self.params.__iter__()

    def __repr__(self):
        return "LaunchParams({:s})".format(repr(self.params))

    def __bool__(self):
        return bool(self.params)

    def __eq__(self, other):
        a, b = self.params, other
        if isinstance(other, LaunchParams):
            b = other.params

        try:
            return len(a) == len(b) and all(a.count(i) == b.count(i) for i in a)
        except:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __getitem__(self, index):
        return self.params[index]

    def __len__(self):
        return len(self.params)

    def union(self, other):
        if self and not other:
            return self

        if other and not self:
            return other

        return LaunchParams(list(itertools.chain(self, other)))

    def product(self, other):
        if self and not other:
            return self

        if other and not self:
            return other

        params = [merge_dict(x, y) for x, y in itertools.product(self, other)]
        return LaunchParams(params)


def newest_output_file(directory):
    """Return path to the newest output file in `directory`.

    Output files obey the convention '*_data-?????.h5'. By newest we
    mean the last one when sorted alphabetically.
    """

    pattern = directory + "/*_data-?????.h5"
    files = sorted(glob.glob(pattern))
    return os.path.abspath(files[-1])


def all_combinations(choices):
    keys = list(choices.keys())
    lp = LaunchParams([{keys[0]: val} for val in choices[keys[0]]])

    for key in keys[1:]:
        lp = lp.product([{key: val} for val in choices[key]])

    return lp


def all_pointwise_combinations(choices):
    n_options = len(choices[next(iter(choices.keys()))])
    assert all(n_options == len(choices[key]) for key in choices)

    lp = LaunchParams()
    for k in range(n_options):
        lp = lp.union(LaunchParams([{key: choices[key][k] for key in choices}]))

    return lp


def pointwise_combinations(choices):
    return all_pointwise_combinations(choices)


def product_of_pointwise_combinations(solvers, pw_choices):
    return solvers.product(all_pointwise_combinations(pw_choices))


def folder_name(scheme):
    return scheme.folder_name()


def build_target(target):
    host = get_host()

    if host == "daint":
        raise Exception("Implement proper build instruction for 'daint'.")

        # cmd = "bash -l -c 'source ~/bin/build_env.sh; make -j12'"
        # subprocess.check_call(cmd, shell=True)

    elif host == "euler":
        raise Exception("Implement proper build instruction for 'euler'.")
        # cmd = "bash -l -c 'source ~/bin/build_env.sh; make -j12'"
        # subprocess.check_call(cmd, shell=True)

    else:
        nproc = multiprocessing.cpu_count()
        cmd = "cd build-release; make -j{} {:s}".format(nproc, target)
        subprocess.check_call(cmd, shell=True)


def build_zisa():
    build_target("zisa")


def build_all():
    build_target("all")
