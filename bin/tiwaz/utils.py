def merge_dict(d1, d2):
    """Merge two dictionaries, i.e. {**d1, **d2} in Python 3.5 onwards."""
    d12 = d1.copy()
    d12.update(d2)
    return d12


def dhhmmss(t):
    days = t.days
    hours, seconds = divmod(t.seconds, 3600)
    minutes, seconds = divmod(seconds, 60)

    return "{:d}-{:02d}:{:02d}:{:02d}".format(days, hours, minutes, seconds)


def hhmm(t):
    days = t.days
    hours, seconds = divmod(t.seconds, 3600)
    minutes, _ = divmod(seconds, 60)
    return "{:02d}:{:02d}".format(24*days + hours, minutes)


def hhmmss(t):
    days = t.days
    hours, seconds = divmod(t.seconds, 3600)
    minutes, seconds = divmod(seconds, 60)
    return "{:02d}:{:02d}:{:02d}".format(24*days + hours, minutes, seconds)


def tail(filename):
    """Perform `tail -1 filename`."""
    with open(filename, 'rb') as f:
        f.seek(-1024, os.SEEK_END)
        return f.readlines()[-1].decode()

def read_txt(filename):
    with open(filename, "r") as f:
        return f.read()

def write_txt(filename, string):
    with open(filename, "w") as f:
        f.write(string)