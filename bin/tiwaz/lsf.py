import subprocess

from .site_details import has_lsf
from .utils import hhmm, read_txt, write_txt


class LSF(object):
    def __init__(self, queue_args):
        assert has_lsf()
        self.queue_args = queue_args

    def submit(self, directory, launch_params, cmd):
        lsf_cmd = self.wrap(launch_params, cmd)

        str_cmd = " ".join([str(c) for c in lsf_cmd])
        with open("runjobs.sh", "a") as f:
            f.write(f"cd {directory}\n")
            f.write(str_cmd + "\n")
            f.write("cd -\n\n")

        print(str_cmd)
        # subprocess.call(cmd, cwd=directory)

    def wrap(self, launch_params, cmd):
        queue_args = self.queue_args

        c = ["bsub", "-J", launch_params.short_id()]

        if hasattr(queue_args, "lsf_args"):
            c += queue_args.lsf_args

        if hasattr(queue_args, "wall_clock"):
            c += ["-W", hhmm(queue_args.wall_clock(launch_params))]

        if hasattr(queue_args, "n_mpi_tasks"):
            # MPI has been requested.
            n_mpi_tasks = queue_args.n_mpi_tasks(launch_params)
            mem = queue_args.memory_per_core(launch_params) * 1e-6

            c += ["-n", str(n_mpi_tasks), "-R", f"rusage[mem={mem:.0f}]"]
            cmd = [f"mpirun -np {n_mpi_tasks} " + " ".join(cmd)]

        elif hasattr(queue_args, "n_omp_threads"):
            # OpenMP has been requested.
            raise Exception("Implement first.")
        else:
            # Okay, must be serial.
            c += ["-n", "1"]

        return c + cmd
