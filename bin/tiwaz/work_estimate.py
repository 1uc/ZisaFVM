import numpy as np
from tiwaz.scheme import read_n_cells


class ZisaWorkEstimate:
    """Estimates the parallelizable work and time required by one CPU. """

    def __init__(self, n0, t0, b0, o0):
        """Estimates work and resources by scaling values from level L.

        Arguments:
            n0: number of cells on level L
            t0: serial runtime on level L
            b0: total memory on level L
            o0: overhead on level L
        """

        self.n0, self.t0 = n0, t0
        self.b0, self.o0 = b0, o0

    def cpu_hours(self, launch_params):
        n_cells = self.n_cells(launch_params)
        return self.t0 * n_cells / self.n0

    def n_cells(self, launch_params):
        return read_n_cells(launch_params["grid"]["file"])

    def work(self, launch_param):
        """Returns an estimate of the total amount of parallelizable work.

        The code is expected to scale up to the point where one core performs
        approximately one unit of work.
        """

        n_cells = self.n_cells(launch_param)
        return n_cells / 1024

    def memory_usage(self, launch_param):
        """Return the required memory.

        This returns `(overhead_per_process, total_memory)` where `overhead_per_process` is the amount
        of 'overhead' each process has and `total_memory` is an estimate of the number of bytes needed
        in a serial run.

        An MPI only run with `n_ranks` ranks requires

           size_per_rank = overhead_per_process + total_memory / n_ranks

        bytes of memory per rank.

        An OMP only run with `n_threads` threads requires

           size_per_thread = (overhead_per_process + total_memory) / n_threads

        bytes of memory per thread.

        An MPI/OMP run with `n_ranks` and `n_threads` per rank, requires the
        same amount of memory as an MPI run with `n_ranks`.
        """
        n_cells = self.n_cells(launch_param)

        overhead_per_process = n_cells / self.n0 * self.o0
        total_memory = n_cells / self.n0 * self.b0

        return overhead_per_process, total_memory


class ZisaFixedMemoryWorkEstimate:
    """Estimates the parallelizable work and time required by one CPU. """

    def __init__(self, n0, t0, b0, o0):
        self.n0, self.t0 = n0, t0
        self.b0, self.o0 = b0, o0

    def cpu_hours(self, launch_params):
        n_cells = self.n_cells(launch_params)
        return self.t0 * n_cells / self.n0

    def n_cells(self, launch_params):
        return read_n_cells(launch_params["grid"]["file"])

    def work(self, launch_param):
        n_cells = self.n_cells(launch_param)
        return n_cells / 1024

    def memory_usage(self, launch_param):
        overhead_per_process = self.o0
        total_memory = self.b0

        return overhead_per_process, total_memory
