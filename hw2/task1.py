import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

xticks = range(10, 31)

times = []
with open("task1.out") as f:
    for i in xticks:
        f.readline()  # Title line
        times.append(float(f.readline()))
        f.readline()  # First element
        f.readline()  # Last element
        f.readline()  # Empty line

with PdfPages("task1.pdf") as pdf:
    plt.plot(xticks, times, "o")
    for x, y in zip(xticks, times):
        if y > 100:
            plt.text(x, y, y)
    plt.xlabel("Array length (2^x)")
    plt.ylabel("Time (milliseconds)")
    plt.xticks(xticks)
    pdf.savefig()
