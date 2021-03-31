from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

begin = 1
end = 10
name = Path(__file__).stem

num_lines_before = 2
num_lines_after = 1

filename_dict = {
    'w/simd': "task2_w_simd.out",
    'wo/simd': "task2_wo_simd.out",
}

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


def read_and_plot(key):
    time = read_times(filename_dict[key])
    plt.plot(xticks, time, "-", label=key)
    for x, y in zip(xticks, time):
        plt.text(
            x, y, f"{y:.3f}",
            horizontalalignment='center',
            fontsize=6,
        )


with PdfPages(f"{name}.pdf") as pdf:
    for key in filename_dict:
        read_and_plot(key)

    plt.title(name)
    plt.xlabel("t")
    plt.ylabel("Time (ms)")
    plt.xticks(xticks)
    plt.legend()

    pdf.savefig()
