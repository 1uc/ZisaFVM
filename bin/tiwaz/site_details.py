import socket
import datetime
import os
import re
import subprocess

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
    def __init__(self, host=None):
        if host is None:
            host = get_host()

        if host == "euler":
            self.cores_per_node = 12
            self.max_nodes = 20
            self.work_per_core = 2.0

        elif host == "daint":
            self.cores_per_node = 12
            self.max_nodes = 160
            self.work_per_core = 2.0

        elif host == "rogui":
            self.cores_per_node = 1
            self.max_nodes = 16
            self.work_per_core = 1.0

        elif host == "liara":
            self.cores_per_node = 1
            self.max_nodes = 1
            self.work_per_core = 1.0

        elif host == "aoifa":
            self.cores_per_node = 1
            self.max_nodes = 12
            self.work_per_core = 1.0

        elif host == "ada":
            nproc = int(subprocess.check_output(["nproc"]))
            self.cores_per_node = 2
            self.max_nodes = nproc // 2
            self.work_per_core = 1.0

        else:
            raise Exception("Unknown host [{:s}]".format(host))

        self.max_cores = self.max_nodes * self.cores_per_node

    def n_tasks(self, work):
        n_proc = int(work / self.work_per_core)

        # Always ask for an entire node.
        cpn = self.cores_per_node
        n_proc = ((n_proc + cpn - 1) // cpn) * cpn
        return max(2, min(n_proc, self.max_cores))

    def memory_per_core(self, overhead_per_process, total_memory, n_tasks):
        return overhead_per_process + total_memory / n_tasks
