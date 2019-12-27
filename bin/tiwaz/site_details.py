import socket
import datetime
import os
import re

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


def zisa_home_directory():
    here = os.path.dirname(os.path.realpath(__file__))
    zisa = os.path.abspath(os.path.join(here, os.pardir, os.pardir))

    return zisa


def get_host():
    """Return the name of the host, strip names like `daint01`, `daint02`."""

    known_hosts = ["daint", "liara", "rogui", "eu-login", "ada", "aoifa"]
    hostname = socket.gethostname()
    match = re.match("^({:s}).*".format("|".join(known_hosts)), hostname)

    if not match:
        raise Exception("Can't deduce hostname from '{:s}'".format(hostname))

    if match.group(1) == "eu-login":
        return "euler"

    return match.group(1)


class MPIHeuristics:
    def __init__(self):
        host = get_host()

        if host == "euler":
            self.cores_per_node = 12
            self.max_nodes = 8
            self.min_block_size = 128

        elif host == "daint":
            self.cores_per_node = 12
            self.max_nodes = 160
            self.min_block_size = 64

        elif host == "rogui":
            self.cores_per_node = 1
            self.max_nodes = 6
            self.min_block_size = 32

        elif host == "liara":
            self.cores_per_node = 1
            self.max_nodes = 1
            self.min_block_size = 32

        elif host == "ada":
            nproc = int(subprocess.check_output(["nproc"]))
            self.cores_per_node = 2
            self.max_nodes = nproc // 2
            self.min_block_size = 32

        else:
            raise Exception("Unknown host [{:s}]".format(get_host()))

        self.max_cores = self.max_nodes * self.cores_per_node

    def n_tasks(self, launch_param):
        n_cells = int(launch_param["nx"]) * int(launch_param["ny"])
        n_proc = int(n_cells / self.min_block_size ** 2)

        # Always ask for an entire node.
        cpn = self.cores_per_node
        n_proc = ((n_proc + cpn - 1) // cpn) * cpn
        return max(1, min(n_proc, self.max_cores))
