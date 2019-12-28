import subprocess
import re


class SLURM(object):
    def __init__(self, queue_args):
        assert has_slurm()
        self.queue_args = queue_args

    def restart(self, launch_param, directory, cmd):
        sbatch_file_pattern = "{:s}/slurm/oneliner-{:s}.sbatch"
        sbatch_file = sbatch_file_pattern.format(directory, get_host())
        cmd = wrap_for_slurm(self.queue_args, sbatch_file, launch_param, cmd)
        self.submit(directory, cmd)

    def submit(self, directory, cmd):
        output = subprocess.check_output(cmd)
        print(output.decode().strip())
        m = re.match(b"Submitted batch job ([0-9]+)$", output)
        pid = m.group(1).decode()

        with open(self.pid_file(directory), "w") as f:
            f.write(pid)

    def pid(self, output_directory):
        with open(self.pid_file(output_directory), "r") as f:
            return f.read()

    def pid_file(self, output_directory):
        return output_directory + "/slurm.pid"

    def is_pid_in_queue(self, pid):
        cmd = ["squeue", "-o", "'%i'"]
        output = subprocess.check_output(cmd)
        output = output.decode().split("\n")
        output = [x[1:-1] for x in output]
        return str(pid) in output

    def wrap(self, sbatch_file, launch_param, cmd):
        queue_args = self.queue_args

        c = ["sbatch"]
        if "job-name" in queue_args:
            c += ["--job-name", queue_args["job-name"]]
        else:
            c += ["--job-name", launch_param["experiment"]]

        if "wall-clock" in queue_args:
            t = queue_args["wall-clock"]
            c += ["--time", dhhmmss(queue_args["wall-clock"])]

        if "n_nodes" in queue_args:
            c += ["-N", str(queue_args["n_nodes"])]

        if "n_tasks" in queue_args:
            c += ["-n", str(queue_args["n_tasks"](launch_param))]

        if "tasks_per_node" in queue_args:
            c += ["--tasks-per-node", str(queue_args["tasks_per_node"])]

        if "slurm_args" in queue_args:
            c += queue_args["slurm_args"]

        return c + [sbatch_file] + cmd
