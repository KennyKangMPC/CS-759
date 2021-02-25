import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

xticks = [2 ** x for x in range(10, 29)]


def read_times(filename):
    times = []
    with open(filename) as f:
        for _ in xticks:
            f.readline()  # Title line
            f.readline()  # Last element
            times.append(float(f.readline()))
            f.readline()  # Empty line
    return times


def read_and_plot(num_threads):
    times = read_times(f"task2-{num_threads}.out")
    plt.plot(xticks, times, "-", label=f"{num_threads} Threads")
    for x, y in zip(xticks, times):
        plt.text(
            x, y, f"{y:.1f}",
            horizontalalignment='center',
            fontsize=6,
        )


with PdfPages("task2.pdf") as pdf:
    read_and_plot(512)
    read_and_plot(1024)

    plt.xlabel("Array Length")
    plt.ylabel("Time (ms)")
    plt.xscale('log', basex=2)
    plt.yscale('log', basey=10)
    plt.xticks(xticks)
    plt.legend()

    pdf.savefig()
