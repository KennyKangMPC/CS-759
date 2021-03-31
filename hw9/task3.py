from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

begin = 1
end = 25
name = Path(__file__).stem

num_lines_before = 1
num_lines_after = 1

xlabel = "n"

xticks = [2 ** x for x in range(begin, end + 1)]


def read_times(filename):
    result = []
    with open(filename) as f:
        for _ in xticks:
            for _ in range(num_lines_before):
                f.readline()
            result.append(float(f.readline()))
            for _ in range(num_lines_after):
                f.readline()
    return result


with PdfPages(f"{name}.pdf") as pdf:
    time = read_times(f"{name}.out")
    plt.plot(xticks, time, "-")
    for x, y in zip(xticks, time):
        plt.text(
            x, y, f"{y:.2f}",
            horizontalalignment='center',
            fontsize=6,
        )

    plt.title(name)
    plt.xlabel(xlabel)
    plt.ylabel("Time (ms)")
    plt.xscale('log', basex=2)
    plt.yscale('log', basey=2)
    plt.xticks(xticks)

    pdf.savefig()
