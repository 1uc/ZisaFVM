import subprocess

from . site_details import has_lsf
from . utils import hhmm, read_txt, write_txt

class LSF(object):
    def __init__(self, queue_args):
        assert(has_lsf())
        self.queue_args = queue_args

    def submit(self, directory, launch_param, cmd):
        lsf_cmd = self.wrap(launch_param, cmd)
        subprocess.call(cmd, cwd=directory)

    def wrap(self, launch_param, cmd):
        queue_args = self.queue_args

        c = ["bsub", "-J", launch_param.short_id()]

        if "lsf_args" in queue_args:
            c += queue_args["lsf_args"]

        if "wall-clock" in queue_args:
            c += ["-W", hhmm(queue_args["wall-clock"])]

        if "n_mpi_tasks" in queue_args:
            # MPI has been requested.
            raise Exception("Implement first.")
        elif "n_omp_threads" in queue_args:
            # OpenMP has been requested.
            raise Exception("Implement first.")
        else:
            # Okay, must be serial.
            c += ["-n", "1"]

        return c + cmd


