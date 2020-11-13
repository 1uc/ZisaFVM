from .site_details import MPIHeuristics
from datetime import timedelta


class QueueArgs:
    def wall_clock(self, launch_param):
        return launch_param["runtime"]


class MPIQueueArgs(QueueArgs):
    def __init__(
        self,
        work_estimate,
        t_min=timedelta(minutes=1),
        t_max=timedelta(days=1),
        heuristics=None,
    ):
        super().__init__()

        self.t_min, self.t_max = t_min, t_max

        self.work_estimate = work_estimate
        self.heuristics = MPIHeuristics() if heuristics is None else heuristics

        self.use_mpi = True

    def n_mpi_tasks(self, launch_param):
        work = self.work_estimate.work(launch_param)
        return self.heuristics.n_tasks(work)

    def wall_clock(self, launch_param):
        n_tasks = self.n_mpi_tasks(launch_param)
        estimate = self.work_estimate.cpu_hours(launch_param) / n_tasks

        if estimate > self.t_max:
            raise RuntimeError(
                f"Warning: estimated runtime exceeds maximum. {estimate}"
            )

        return min(max(estimate, self.t_min), self.t_max)

    def memory_per_core(self, launch_param):
        memory_usage = self.work_estimate.memory_usage(launch_param)
        overhead_per_process, total_memory = memory_usage

        n_tasks = self.n_mpi_tasks(launch_param)
        sharp_requirement = overhead_per_process + total_memory / n_tasks
        return sharp_requirement


class FixedMPIQueueArgs(QueueArgs):
    def __init__(self, mem_per_core, wall_clock):
        self._mem_per_core = mem_per_core
        self._wall_clock = wall_clock

        self.use_mpi = True

    def n_mpi_tasks(self, launch_params):
        return launch_params["experiment"]["n_proc"]

    def wall_clock(self, launch_params):
        return self._wall_clock

    def memory_per_core(self, launch_params):
        return self._mem_per_core
