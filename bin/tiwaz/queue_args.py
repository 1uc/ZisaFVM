from .site_details import MPIHeuristics
from datetime import timedelta


class QueueArgs:
    def wall_clock(self, launch_param):
        return launch_param["runtime"]


class MPIQueueArgs(QueueArgs):
    def __init__(
        self, work_estimate, t_min=timedelta(minutes=1), t_max=timedelta(days=1)
    ):
        super().__init__()

        self.t_min, self.t_max = t_min, t_max

        self.work_estimate = work_estimate
        self.heuristics = MPIHeuristics()
        self.use_mpi = True

    def n_mpi_tasks(self, launch_param):
        work = self.work_estimate.work(launch_param)
        return self.heuristics.n_tasks(work)

    def wall_clock(self, launch_param):
        n_tasks = self.n_mpi_tasks(launch_param)
        estimate = self.work_estimate.cpu_hours(launch_param) / n_tasks

        if estimate > self.t_max:
            raise RuntimeError("Warning: estimated runtime exceeds maximum.")

        return min(max(estimate, self.t_min), self.t_max)

    def memory_per_core(self, launch_param):
        memory_usage = self.work_estimate.memory_usage(launch_param)
        n_tasks = self.n_mpi_tasks(launch_param)

        sharp_requirement = self.heuristics.memory_per_core(*memory_usage, n_tasks)
        return 1.5 * sharp_requirement
