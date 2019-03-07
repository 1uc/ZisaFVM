import subprocess

from . site_details import has_no_queue

class NoQueue:
    def __init__(self, queue_args):
        assert(has_no_queue())
        self.queue_args = queue_args

    def submit(self, directory, cmd):
        subprocess.call(cmd, cwd=directory)
