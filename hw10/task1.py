import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

begin = 1
end = 10

num_lines_before = 1
num_lines_after = 0

xticks = ["optimize1", "optimize2", "optimize3", "optimize4", "optimize5", "opt1-SIMD"]


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


for name in ["task11", "task12", "task13", "task14"]:
    with PdfPages(f"{name}.pdf") as pdf:
        time = read_times(f"{name}.out")
        plt.plot(xticks, time, "-")
        for x, y in zip(xticks, time):
            plt.text(
                x, y, f"{y:.4f}",
                horizontalalignment='center',
            )

        plt.title(name)
        plt.ylabel("Time (ms)")
        plt.xticks(xticks)

        pdf.savefig()

        plt.clf()
