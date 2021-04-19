import socket
import datetime
import os
import re
import subprocess
import toml

import scibs

hosts_with_slurm = ["daint"]
hosts_with_lsf = ["euler"]
hosts_without_queue = ["rogui", "liara", "ada", "aoifa"]


def has_slurm():
    return get_host() in hosts_with_slurm


def has_lsf():
    return get_host() in hosts_with_lsf


def has_no_queue():
    return get_host() in hosts_without_queue


def todays_scratch_directory():
    scratch = "/".join(["${SCRATCH}", datetime.date.today().isoformat()])
    scratch = os.path.expandvars(scratch)
    latest = os.path.expandvars("${SCRATCH}" + "/latest")

    if os.path.islink(latest):
        os.unlink(latest)
    os.symlink(scratch, latest)

    return scratch


def zisa_config():
    here = os.path.dirname(os.path.realpath(__file__))
    tiwaz = os.path.abspath(os.path.join(here, os.pardir, os.pardir))

    config_file = os.path.join(tiwaz, "config.toml")
    if os.path.exists(config_file):
        return toml.load(config_file)

    return None


def zisa_install_directory():
    config = zisa_config()
    if "zisa" in config and "install_directory" in config["zisa"]:
        return os.path.expandvars(
            os.path.expanduser(config["zisa"]["install_directory"])
        )

    else:
        here = os.path.dirname(os.path.realpath(__file__))
        return os.path.abspath(os.path.join(here, os.pardir, os.pardir, "zisa"))


def zisa_home_directory():
    config = zisa_config()
    if "zisa" in config and "home_directory" in config["zisa"]:
        return os.path.expanduser(os.path.expandvars(config["zisa"]["home_directory"]))

    raise Exception("Unknown Zisa root directory.")


def zisa_build_directory():
    config = zisa_config()
    if "zisa" in config and "build_directory" in config["zisa"]:
        return os.path.expanduser(os.path.expandvars(config["zisa"]["build_directory"]))

    raise Exception("Unknown Zisa root directory.")


def zisa_grid_repository():
    default_grid_repository = os.path.join(zisa_home_directory(), "grids")

    try:
        host = get_host()
        if host in ["euler", "leonhard"]:
            return os.path.expandvars("${WORK}/grids")
        else:
            scratch = os.path.expandvars("${SCRATCH}")
            if scratch != "":
                return os.path.join(scratch, "grids")
            else:
                return default_grid_repository

    except UnknownHostError:
        return default_grid_repository


class UnknownHostError(Exception):
    """Indicates that the 'generic hostname' could not be determined."""

    pass


def get_host():
    """Return the name of the host, strip names like `daint01`, `daint02`."""

    known_hosts = ["daint", "liara", "rogui", "eu-login", "ada", "aoifa"]
    hostname = socket.gethostname()
    match = re.match("^({:s}).*".format("|".join(known_hosts)), hostname)

    if not match:
        raise UnknownHostError("Can't deduce hostname from '{:s}'".format(hostname))

    if match.group(1) == "eu-login":
        return "euler"

    return match.group(1)


class MPIHeuristics:
    def __init__(
        self, host=None, cores_per_node=None, max_nodes=None, work_per_core=None
    ):

        if all(x is not None for x in [cores_per_node, max_nodes, work_per_core]):
            self._init_by_hand(cores_per_node, max_nodes, work_per_core)

        else:
            if host is None:
                host = get_host()

            self._init_by_host(host)

        self.max_cores = self.max_nodes * self.cores_per_node

    def _init_by_hand(self, cores_per_node, max_nodes, work_per_core):
        self.cores_per_node = cores_per_node
        self.max_nodes = max_nodes
        self.work_per_core = work_per_core

    def _init_by_host(self, host):
        if host == "euler":
            self.max_cores_per_node = 128
            self.cores_per_node = 12
            self.max_nodes = 20
            self.work_per_core = 2.0

        elif host == "daint":
            self.cores_per_node = 12
            self.cores_per_node = 12
            self.max_nodes = 160
            self.work_per_core = 2.0

        elif host == "rogui":
            self.max_cores_per_node = 16
            self.cores_per_node = 1
            self.max_nodes = 16
            self.work_per_core = 1.0

        elif host == "liara":
            self.max_cores_per_node = 2
            self.cores_per_node = 1
            self.max_nodes = 2
            self.work_per_core = 1.0

        elif host == "aoifa":
            self.max_cores_per_node = 12
            self.cores_per_node = 1
            self.max_nodes = 12
            self.work_per_core = 1.0

        elif host == "ada":
            nproc = int(subprocess.check_output(["nproc"]))
            self.max_cores_per_node = nproc
            self.cores_per_node = 2
            self.max_nodes = nproc // 2
            self.work_per_core = 1.0

        else:
            raise Exception("Unknown host [{:s}]".format(host))

    def n_tasks(self, work):
        n_proc = int(work / self.work_per_core)

        # Always ask for an entire node.
        cpn = self.cores_per_node
        n_proc = ((n_proc + cpn - 1) // cpn) * cpn
        return max(2, min(n_proc, self.max_cores))


def max_cores_per_node():
    heuristics = MPIHeuristics(host)
    return heuristics.max_cores_per_node


def make_batch_system():
    host = get_host()
    if host == "euler":
        return scibs.EulerLSF()

    elif has_slurm():
        return scibs.SLUM()

    elif has_lsf():
        return scibs.LSF()

    else:
        return scibs.LocalBS()
