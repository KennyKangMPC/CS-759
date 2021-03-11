import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

xticks = [2 ** x for x in range(2, 16)]


def read_times(filename):
    result = []
    with open(filename) as f:
        for _ in xticks:
            f.readline()  # Title line
            result.append(float(f.readline()))
            f.readline()  # Empty line
    return result


times = {
    'with tensor': read_times(f"task1-w-tensor.out"),
    'without tensor': read_times(f"task1-wo-tensor.out")
}

diff = {x: t1 - t2 for x, t1, t2 in zip(xticks, *times.values())}


def get_vertical_alignment(x, key):
    if key == list(times)[0]:
        return 'top' if diff[x] < 0 else 'bottom'
    else:
        return 'top' if diff[x] > 0 else 'bottom'


def read_and_plot(key):
    plt.plot(xticks, times[key], "-", label=key)
    for x, y in zip(xticks, times[key]):
        plt.text(
            x, y, f"{y:.1f}",
            horizontalalignment='center',
            verticalalignment=get_vertical_alignment(x, key),
            fontsize=6,
        )


with PdfPages("task1.pdf") as pdf:
    for key in times:
        read_and_plot(key)

    plt.xlabel("n")
    plt.ylabel("Time (ms)")
    plt.xscale('log', basex=2)
    plt.yscale('log', basey=10)
    plt.xticks(xticks)
    plt.legend()

    pdf.savefig()
