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

        is_mpi = hasattr(queue_args, "n_mpi_tasks")
        is_omp = hasattr(queue_args, "n_omp_threads")

        if is_mpi and not is_omp:
            # MPI has been requested.
            n_mpi_tasks = queue_args.n_mpi_tasks(launch_params)
            mem = queue_args.memory_per_core(launch_params) * 1e-6
            # ptile = min(128, n_mpi_tasks)  # FIXME this is site specific

            c += [
                "-n",
                str(n_mpi_tasks),
                "-R",
                f"rusage[mem={mem:.0f}]",
                # "-R",
                # f"span[ptile={ptile}]",
            ]
            cmd = [f"mpirun -np {n_mpi_tasks} " + " ".join(cmd)]

        elif is_mpi and is_omp:
            # OpenMP has been requested.
            n_mpi_tasks = queue_args.n_mpi_tasks(launch_params)
            n_omp_threads = queue_args.n_omp_threads(launch_params)
            n_cores = n_mpi_tasks * n_omp_threads
            mem = queue_args.memory_per_core(launch_params) * 1e-6
            ptile = min(128, n_cores)  # FIXME site-specific

            c += [
                "-n",
                str(n_cores),
                "-R",
                f"rusage[mem={mem:.0f}]",
                "-R",
                f"span[ptile={ptile}]",
            ]
            cmd = [
                " ".join(
                    [
                        f"'export OMP_NUM_THREADS={n_omp_threads};",
                        "unset LSF_AFFINITY_HOSTFILE;",
                        f"mpirun -np {n_mpi_tasks} --map-by node:PE={n_omp_threads}",
                    ]
                    + cmd
                    + ["'"]
                )
            ]

        elif not is_mpi and is_omp:
            # OpenMP has been requested.
            raise Exception("Implement first.")

        else:
            # Okay, must be serial.
            c += ["-n", "1"]

        return c + cmd
