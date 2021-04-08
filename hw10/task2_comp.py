from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

name = Path(__file__).stem

num_lines_before = 2
num_lines_after = 1

filename_dict = {
    'OpenMP & MPI': "task2_2.out",
    'Pure OpenMP': "task2_pure_omp_2.out",
}

xticks = [2 ** x for x in range(1, 27)]


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


def read_and_plot(key):
    time = read_times(filename_dict[key])
    plt.plot(xticks, time, "-", label=key)
#    for x, y in zip(xticks, time):
#        plt.text(
#            x, y, f"{y:.f}",
#            horizontalalignment='center',
#            fontsize=6,
#        )


with PdfPages(f"{name}.pdf") as pdf:
    for key in filename_dict:
        read_and_plot(key)

    plt.title("TASK2D2 t=9")
    plt.xlabel("n")
    plt.ylabel("Time (ms)")
    plt.xscale('log', basex=2)
    plt.yscale('log', basey=2)
    plt.xticks(xticks)
    plt.legend()

    pdf.savefig()
