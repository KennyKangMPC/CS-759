import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

xticks = [2 ** x for x in range(10, 30)]


def read_times(filename):
    times = []
    with open(filename) as f:
        for _ in xticks:
            f.readline()  # Title line
            times.append(float(f.readline()))
            f.readline()  # First element
            f.readline()  # Last element
            f.readline()  # Empty line
    return times


def read_and_plot(num_threads):
    times = read_times(f"task3-{num_threads}.out")
    plt.plot(xticks, times, ".", label=f"{num_threads} Threads")
    for x, y in zip(xticks, times):
        plt.text(
            x, y, f"{y:.4g}",
            horizontalalignment='center',
            fontsize=6,
        )


with PdfPages("task3.pdf") as pdf:
    read_and_plot(16)
    read_and_plot(512)

    plt.xlabel("Array Size (2^x)")
    plt.ylabel("Time (milliseconds)")
    plt.xscale('log', basex=2)
    plt.yscale('log', basey=10)
    plt.title("Scaling Analysis Plot")
    plt.xticks(xticks)
    plt.legend()

    pdf.savefig()
