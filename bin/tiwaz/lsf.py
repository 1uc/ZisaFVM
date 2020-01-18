import subprocess

from .site_details import has_lsf
from .utils import hhmm, read_txt, write_txt


class LSF(object):
    def __init__(self, queue_args):
        assert has_lsf()
        self.queue_args = queue_args

    def submit(self, directory, launch_param, cmd):
        lsf_cmd = self.wrap(launch_param, cmd)
        print(" ".join(lsf_cmd))
        # subprocess.call(cmd, cwd=directory)

    def wrap(self, launch_param, cmd):
        queue_args = self.queue_args

        c = ["bsub", "-J", launch_param.short_id()]

        if hasattr(queue_args, "lsf_args"):
            c += queue_args["lsf_args"]

        if hasattr(queue_args, "wall_clock"):
            c += ["-W", hhmm(queue_args.wall_clock(launch_param))]

        if hasattr(queue_args, "n_mpi_tasks"):
            # MPI has been requested.
            c += ["-n", str(queue_args.n_mpi_tasks(launch_param))]
            cmd = ["mpirun " + " ".join(cmd)]

        elif hasattr(queue_args, "n_omp_threads"):
            # OpenMP has been requested.
            raise Exception("Implement first.")
        else:
            # Okay, must be serial.
            c += ["-n", "1"]

        return c + cmd
