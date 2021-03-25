from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

begin = 1
end = 20
name = Path(__file__).stem

num_lines_before = 3
num_lines_after = 1

xticks = range(begin, end + 1)


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
            x, y, f"{y:.1f}",
            horizontalalignment='center',
            fontsize=6,
        )

    plt.title("Plot of task3_t")
    plt.xlabel("t")
    plt.ylabel("Time (ms)")
    plt.xticks(xticks)

    pdf.savefig()
