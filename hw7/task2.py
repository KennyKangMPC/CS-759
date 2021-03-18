import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

xticks = [2 ** x for x in range(5, 25)]


def read_times(filename):
    result = []
    with open(filename) as f:
        for _ in xticks:
            f.readline()  # Title line
            f.readline()  # Result line
            f.readline()  # Result line
            result.append(float(f.readline()))
            f.readline()  # Empty line
    return result


with PdfPages("task2.pdf") as pdf:
    time = read_times(f"task2.out")
    plt.plot(xticks, time, "-")
    for x, y in zip(xticks, time):
        plt.text(
            x, y, f"{y:.2f}",
            horizontalalignment='center',
            fontsize=6,
        )

    plt.xlabel("n")
    plt.ylabel("Time (ms)")
    plt.xscale('log', basex=2)
    plt.yscale('log', basey=2)
    plt.title("Scaling Analysis Plot")
    plt.xticks(xticks)

    pdf.savefig()
