import shutil
import glob

from . site_details import *
from . queue import make_queue
from . launch_params import folder_name
from . post_process import find_data_files, find_grid


class LaunchJob(object):
    """Base class for launching jobs."""

    def __init__(self, base_directory, queue_args = None):
        """Create a `LaunchJob`.

        :base_directory: store the output in subfolders of `base_directory`.
        :queue_args: arguments to be passed to the queue.
        """

        self.queue = make_queue(queue_args if queue_args else [])
        self.base_directory = base_directory

    def launch_command(self, launch_param, directory):
        return self.basic_launch_command(launch_param, directory)

    def basic_launch_command(self, launch_param, directory):
        binary_name = directory + "/zisa"
        config_name = directory + "/config.json"
        cmd = [binary_name, "run", "--config", config_name]

        return cmd

    def output_directory(self, launch_param):
        directory = "/".join([self.base_directory, folder_name(launch_param)])
        return os.path.expandvars(directory)

    def log_file(self, directory):
        return directory + "/zisa.log"


class LaunchNewJob(LaunchJob):
    """Takes care of everything required to launch a new simulation."""

    def __init__(self, zisa_home, base_directory, force, queue_args = None):
        super().__init__(base_directory, queue_args)

        self.zisa_home = zisa_home
        self.force = force

    def __call__(self, launch_param):
        directory = self.output_directory(launch_param)
        self.preparatory_work(launch_param)
        cmd = self.launch_command(launch_param, directory)

        self.queue.submit(directory, launch_param, cmd)

    def preparatory_work(self, launch_param):
        directory = self.make_clean_directory(launch_param)
        self.copy_essentials(directory, launch_param)

    def make_clean_directory(self, launch_param):
        directory = self.output_directory(launch_param)

        if self.force and os.path.isdir(directory):
            shutil.rmtree(directory)

        os.makedirs(directory)
        return directory

    def copy_essentials(self, directory, launch_param):
        zisa_home = self.zisa_home

        self.copy_binaries(zisa_home, directory)
        self.write_config(directory, launch_param)
        self.copy_grid(zisa_home, directory, launch_param)
        self.copy_shaders(zisa_home, directory)
        self.write_gitinfo(zisa_home, directory)

    def copy_binaries(self, zisa_home, directory):
        self.copy(zisa_home + "/build-release", directory, ["zisa"])

    def copy_grid(self, zisa_home, directory, launch_params):
        grid_paths = [launch_params.grid_filename()]

        if "reference" in launch_params:
            grid_paths += launch_params["reference"]["coarse_grids"]

        self.copy(zisa_home, directory, grid_paths)

    def copy_shaders(self, zisa_home, directory):
        shaders = glob.glob("shaders/*")
        self.copy(zisa_home, directory, shaders)

    def write_config(self, directory, launch_param):
        launch_param.save(directory + "/config.json")

    def write_gitinfo(self, zisa_home, directory):
        print("warning: remember to implement git info.")
        # self.copy(zisa_home + "/build", directory, ["git.diff", "git.hash", "git.status"])

    def copy(self, src, dest, files):
        for f in files:
            filename = os.path.join(dest, f)
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            shutil.copy(os.path.join(src, f), filename)

class RestartJob(LaunchJob):
    """Takes care of everything required to launch a new simulation."""

    def __init__(self, base_directory, queue_args = None):
        super().__init__(base_directory, queue_args)

    def __call__(self, launch_param):
        directory = self.output_directory(launch_param)
        cmd = self.launch_command(launch_param, directory)

        self.write_config(directory, launch_param)
        self.queue.submit(directory, launch_param, cmd)

    def write_config(self, directory, launch_param):
        launch_param.save(directory + "/config.json")


def launch_all(launch_params, force, queue_args=None):
    launch_job = LaunchNewJob(zisa_home_directory(),
                              todays_scratch_directory(),
                              force,
                              queue_args)

    for param in launch_params:
        launch_job(param)

def restart_all(launch_params, restart_index, queue_args=None):
    cwd = os.getcwd()
    restart_job = RestartJob(cwd, queue_args)

    for param in launch_params:
        dir = folder_name(param)
        data_files = find_data_files(dir)
        data_file = os.path.relpath(data_files[restart_index], dir)
        param["restart"] = {"file": data_file}
        param["grid"]["file"] = os.path.relpath(find_grid(dir), dir)

        restart_job(param)
