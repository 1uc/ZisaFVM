import subprocess

from . site_details import has_lsf
from . utils import hhmm, read_txt, write_txt

class LSF(object):
    def __init__(self, queue_args):
        assert(has_lsf())
        self.queue_args = queue_args

    def submit(self, directory, launch_param, cmd):
        self.generate_bsub_file(directory, launch_param)
        lsf_cmd = self.wrap(launch_param, cmd)
        subprocess.call(cmd, cwd=directory)

    def generate_bsub_file(self, directory, launch_param):
        bsub = ""

        bsub = "\n".join([
            self.job_name_line(launch_param),
            self.wall_clock_line()
        ])

    def job_name_line(self, launch_param):
        queue_args = self.queue_args

        if "job-name" in queue_args:
            job_name = queue_args["job-name"]
        else:
            job_name = launch_param["experiment"].short_id()

        return "#BSUB -J {}".format(job_name)

    def wall_clock_line(self):
        wall_clock = hhmm(self.queue_args["wall-clock"])
        return "#BSUB -W {}".format(wall_clock)


    def wrap(self, launch_param, cmd):
        queue_args = self.queue_args

        c = ["bsub"]

        if "lsf_args" in queue_args:
            c += queue_args["lsf_args"]

        if "n_tasks" in queue_args:
            n_tasks = queue_args["n_tasks"](launch_param)
            c += ["-n", str(n_tasks)]
            cmd = ["mpirun", "-n", str(n_tasks)] + cmd

        return c + cmd


