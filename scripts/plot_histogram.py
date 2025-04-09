import struct
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import sys


def read_doubles_from_binary_file(filename):
    doubles = []
    with open(filename, "rb") as f:
        while True:
            data = f.read(8)  # 8 bytes for a double
            if not data:
                break
            doubles.append(struct.unpack("d", data))  # '<d' for little-endian double
    return doubles


def plot_histogram(doubles, xlabel, save_name, xmax, ymax):
    print("we have " + str(len(doubles)) + " data points")
    num_bins = 1000
    counts, bin_edges = np.histogram(doubles, bins=num_bins)
    # Compute bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_widths = bin_edges[1:] - bin_edges[:-1]
    # normalize the data
    normalized = [float(counts[i]) / bin_widths[i] for i in range(len(counts))]
    normalized = counts / bin_widths
    normalized /= len(doubles)
    # Plot histogram as line graph using matplotlib
    plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.xlim(0, xmax)
    plt.ylim(0, ymax)
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
    plt.gcf().canvas.get_default_filename = lambda: save_name
    plt.plot(bin_centers, normalized)
    plt.savefig("full_diagnostics/" + save_name)
    plt.show()


# replace 'file.bin' with your binary file
if len(sys.argv) != 6:
    print(str(len(sys.argv)) + "parameters, expected 4, usage:")
    print(
        "python3 plot_histogram [data_loc] [x_axis_label] [save_file_name] [x_axis_max] [y_axis_max]"
    )
    sys.exit()
doubles = read_doubles_from_binary_file(sys.argv[1])
plot_histogram(
    doubles, sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5])
)
print("done!")
