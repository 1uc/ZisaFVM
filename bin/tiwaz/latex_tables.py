import itertools
import numpy as np

from . post_process import extract_solver_data

class LatexConvergenceTable(object):
    """Format convergence data as LaTeX `tabular`."""

    def __init__(self, x_labels, data):
        self.table = "".join([self.header(x_labels),
                              self.content(data),
                              self.footer()])

    def write(self, filename):
        with open(filename, "w+") as f:
            f.write(self.table)

    def header(self, x_labels):
        n = len(x_labels)

        header = "".join(
            ["\\begin{tabular}{l\n",
            n*"                S[table-format=3.2e2]\n",
            "}\n",
            "\\toprule\n",
            "$\Delta x$ & ",
            " & ".join("\multicolumn{1}{S[table-format=3.2e2]}{" + l + "}" for l in x_labels),
            " \\\\\n",
            "\\midrule\n"])
        return header

    def error_subheader(self):
        return "\n\\multicolumn{1}{c}{err}"

    def rate_subheader(self):
        return "\n\\multicolumn{1}{c}{rate}"

    def footer(self):
        return "".join(["\\bottomrule\n", "\\end{tabular}"])

    def content(self, data):
        content = ""
        for k, row in enumerate(data):
            content += " & ".join(row) + "  \\\\\n"

            if k % 2 == 1 and k != len(data)-1:
                content += "\midrule\n"

        return content


def convergence_rate(res, err):
    log_err = np.log(err)
    log_res = np.log(res)
    return (log_err[1:] - log_err[:-1])/(log_res[1:] - log_res[:-1])


def format_table_contents(solver_keys, resolutions, l1_errors, rates):
    return format_x_err_y_solver_rate(solver_keys, resolutions, l1_errors, rates)

def format_x_err_y_solver_rate(solver_keys, resolutions, l1_errors, rates):
    max_str_len = 32
    dtype = "<U" + str(max_str_len)
    table = np.empty((2*len(solver_keys), 1 + len(resolutions)), dtype=dtype)

    width = max(len(key) for key in solver_keys) + 2
    assert(width <= max_str_len)

    table[:-1:2,0] = [key.ljust(width) for key in solver_keys]
    table[1::2,0] = "rate"

    for k, l1_err in enumerate(l1_errors):
        table[2*k,1:] = ["{:8.2e}".format(err) for err in l1_err]

    table[1::2,1] = "\\multicolumn{1}{c}{--}"
    for k, rate in enumerate(rates):
        table[2*k+1,2:] = ["{:8.2f}".format(r) for r in rate]

    return table

def format_x_err_rate_y_solver(solver_keys, resolutions, l1_errors, rates):
    max_str_len = 32
    dtype = "<U" + str(max_str_len)
    table = np.empty((len(solver_keys), 1 + 2*len(resolutions)), dtype=dtype)

    width = max(len(key) for key in solver_keys) + 2
    assert(width <= max_str_len)

    table[:,0] = [key.ljust(width) for key in solver_keys]

    for k, l1_err in enumerate(l1_errors):
        table[k,1:-1:2] = ["{:8.2e}".format(err) for err in l1_err]

    table[:,2] = "--"
    for k, rate in enumerate(rates):
        table[k,4::2] = ["{:8.2f}".format(r) for r in rate]

    return table


def write_convergence_table(results, columns, labels, filename):
    """Write a latex table of convergence rates to disk."""

    for key in ["l1_error", "l1_eq_error"]:
        all_errors = []
        all_rates = []
        solver_keys = [labels(col) for col in columns]

        for col in columns:
            result = extract_solver_data(col, results)
            resolutions = np.array([r["dx_max"] for r in result])
            l1_err = [r[key] for r in result]

            rate = convergence_rate(resolutions, l1_err)

            all_errors.append(l1_err)
            all_rates.append(rate)

        l1_errors = np.array(all_errors)
        rates = np.array(all_rates)

        data = format_table_contents(solver_keys, resolutions, l1_errors, rates)
        x_labels = ["{:8.2e}".format(res) for res in resolutions]

        table = LatexConvergenceTable(x_labels, data)
        table.write(filename + "_" + key + ".tex")
