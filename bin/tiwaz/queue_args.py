from .site_details import MPIHeuristics


class QueueArgs:
    def wall_clock(self, launch_param):
        return launch_param["runtime"]


class MPIQueueArgs(QueueArgs):
    def __init__(self, work_estimate):
        super().__init__()

        self.work_estimate = work_estimate
        self.heuristics = MPIHeuristics()
        self.use_mpi = True

    def n_mpi_tasks(self, launch_param):
        work = self.work_estimate.work(launch_param)
        return self.heuristics.n_tasks(work)

    def wall_clock(self, launch_param):
        n_tasks = self.n_mpi_tasks(launch_param)
        return self.work_estimate.cpu_hours(launch_param) / n_tasks
