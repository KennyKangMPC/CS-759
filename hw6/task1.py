import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

xticks = [2 ** x for x in range(5, 15)]


def read_times(filename):
    result = []
    with open(filename) as f:
        for _ in xticks:
            f.readline()  # Title line
            result.append(float(f.readline()))
            f.readline()  # Empty line
    return result


with PdfPages("task1.pdf") as pdf:
    times = read_times(f"task1.out")
    plt.plot(xticks, times, "-")
    for x, y in zip(xticks, times):
        plt.text(
            x, y, f"{y:.3f}",
            horizontalalignment='center',
        )

    plt.xlabel("n")
    plt.ylabel("Time (ms)")
    plt.xscale('log', basex=2)
    #plt.yscale('log', basey=10)
    plt.title("Scaling Analysis Plot")
    plt.xticks(xticks)

    pdf.savefig()
