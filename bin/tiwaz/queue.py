from .slurm import SLURM
from .no_queue import NoQueue
from .site_details import has_lsf, has_slurm, has_no_queue


def make_queue(queue_args):
    if has_slurm():
        return SLURM(queue_args)

    if has_lsf():
        return LSF(queue_args)

    if has_no_queue():
        return NoQueue(queue_args)

    raise Exception("Can't determine the queue to use.")
