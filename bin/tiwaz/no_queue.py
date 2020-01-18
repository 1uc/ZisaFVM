import subprocess
import os

from .site_details import has_no_queue


class NoQueue:
    def __init__(self, queue_args):
        assert has_no_queue()
        self.queue_args = queue_args

    def submit(self, directory, launch_param, cmd):
        wrapped_command = self.wrap(directory, launch_param, cmd)

        # raw_cmd = " ".join([str(c) for c in wrapped_command])
        # with open("runjobs.sh", "a") as f:
        #     f.write(f"cd {directory} && {raw_cmd} ; cd -")

        subprocess.run(wrapped_command, cwd=directory, check=True)

    def wrap(self, directory, launch_param, cmd):
        if getattr(self.queue_args, "use_mpi", False):
            n_pe = self.queue_args.n_mpi_tasks(launch_param)
            return ["mpirun", "-np", str(n_pe)] + cmd

        else:
            return cmd
