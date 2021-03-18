import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

xticks = [2 ** x for x in range(10, 31)]


def read_times(filename):
    result = []
    with open(filename) as f:
        for _ in xticks:
            f.readline()  # Title line
            f.readline()  # Result line
            result.append(float(f.readline()))
            f.readline()  # Empty line
    return result


times = {
    'CUB': read_times(f"task1_cub.out"),
    'Thrust': read_times(f"task1_thrust.out"),
    'hw5 GPU-1024 threads': read_times(f"../hw5/task1-1024.out"),
    'hw5 GPU-512 threads': read_times(f"../hw5/task1-512.out")
}


def read_and_plot(key):
    plt.plot(xticks, times[key], "-", label=key)
    for x, y in zip(xticks, times[key]):
        plt.text(
            x, y, f"{y:.1f}",
            horizontalalignment='center',
            fontsize=6,
        )


with PdfPages("task1.pdf") as pdf:
    for key in times:
        read_and_plot(key)

    plt.xlabel("n")
    plt.ylabel("Time (ms)")
    plt.xscale('log', basex=2)
    plt.yscale('log', basey=2)
    plt.xticks(xticks)
    plt.title("Scaling Analysis Plot")
    plt.legend()
    pdf.savefig()
