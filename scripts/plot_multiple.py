import struct
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import sys
from collections import defaultdict


def read_labeled_doubles_from_binary_file(filename, labelints):
    record_format = f"{labelints}id"

    with open(filename, "rb") as f:
        data = f.read()

    doubles = defaultdict(list)
    for record in struct.iter_unpack(record_format, data):
        ints = record[:-1]  # First N values are label integers
        value = record[-1]  # Last value is the double
        doubles[ints].append(value)

    return doubles


def plot_histogram(doubles, key, xmax):
    counts, bin_edges = np.histogram(doubles, bins=50, range=(0, xmax))
    # Compute bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_widths = bin_edges[1:] - bin_edges[:-1]
    # normalize the data
    normalized = [float(counts[i]) / bin_widths[i] for i in range(len(counts))]
    normalized = counts / bin_widths
    normalized /= len(doubles)
    # Plot histogram as line graph using matplotlib
    str_key = "-".join(str(num) for num in key)
    plt.plot(bin_centers, normalized, label=str_key)


# replace 'file.bin' with your binary file
if len(sys.argv) != 6:
    print(str(len(sys.argv)) + "parameters, expected 4, usage:")
    print(
        "python3 plot_histogram [data_loc] [x_axis_label] [save_file_name] [x_axis_max] [y_axis_max]"
    )
    sys.exit()
doubles = read_labeled_doubles_from_binary_file(sys.argv[1], 2)
top_items = sorted(doubles.items(), key=lambda item: len(item[1]), reverse=True)[:10]
plt.xlabel(sys.argv[2])
plt.ylabel("Frequency")
plt.xlim(0, float(sys.argv[4]))
plt.ylim(0, float(sys.argv[5]))
current_date = dt.datetime.now().strftime("%Y-%m-%d")
plt.text(
    0.98,
    0.98,
    f"{current_date}",
    transform=plt.gca().transAxes,
    fontsize=10,
    verticalalignment="top",
    horizontalalignment="right",
)
plt.gcf().canvas.get_default_filename = lambda: sys.argv[3]
for item in top_items:
    plot_histogram(item[1], item[0], float(sys.argv[4]))
plt.legend()
font = {"family": "normal", "weight": "bold", "size": 22}

plt.rc("font", **font)
plt.savefig("../plots/" + sys.argv[3])
plt.show()
print("done!")
