import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

xticks = [2 ** x for x in range(10, 31)]


def read_times(filename):
    result = []
    with open(filename) as f:
        for _ in xticks:
            f.readline()  # Title line
            f.readline()  # Sum
            result.append(float(f.readline()))
            f.readline()  # Empty line
    return result


times = {
    512: read_times(f"task1-512.out"),
    1024: read_times(f"task1-1024.out")
}

diff = {x: t1 - t2 for x, t1, t2 in zip(xticks, times[512], times[1024])}


def get_vertical_alignment(x, num_threads):
    if num_threads == 512:
        return 'top' if diff[x] < 0 else 'bottom'
    else:
        return 'top' if diff[x] > 0 else 'bottom'


def read_and_plot(num_threads):
    plt.plot(xticks, times[num_threads], "-", label=f"{num_threads} threads per block")
    for x, y in zip(xticks, times[num_threads]):
        plt.text(
            x, y, f"{y:.1f}",
            horizontalalignment='center',
            verticalalignment=get_vertical_alignment(x, num_threads),
            fontsize=6,
        )


with PdfPages("task1.pdf") as pdf:
    read_and_plot(512)
    read_and_plot(1024)

    plt.xlabel("N")
    plt.ylabel("Time (ms)")
    plt.xscale('log', basex=2)
    plt.yscale('log', basey=10)
    plt.xticks(xticks)
    plt.legend()

    pdf.savefig()
