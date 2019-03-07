class LSF(object):
    def __init__(self):
        self.queue_args = None
        raise Exception("Implement remaining functionality.")

    def wrap(self, launch_param, cmd):
        queue_args = self.queue_args

        c = ["bsub"]

        if "job-name" in queue_args:
            c += ["-J", queue_args["job-name"]]
        else:
            c += ["-J", launch_param["experiment"]]

        if "wall-clock" in queue_args:
            c += ["-W", hhmm(queue_args["wall-clock"])]

        if "n_nodes" in queue_args:
            raise Exception("This won't work...")

        if "tasks_per_node" in queue_args:
            raise Exception("This won't work...")

        if "lsf_args" in queue_args:
            c += queue_args["lsf_args"]

        if "n_tasks" in queue_args:
            n_tasks = queue_args["n_tasks"](launch_param)
            c += ["-n", str(n_tasks)]
            cmd = ["mpirun", "-n", str(n_tasks)] + cmd

        return c + cmd


