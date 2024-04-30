import struct
import matplotlib.pyplot as plt
import numpy as np


def read_doubles_from_binary_file(filename):
    doubles = []
    with open(filename, "rb") as f:
        while True:
            data = f.read(8)  # 8 bytes for a double
            if not data:
                break
            doubles.append(struct.unpack("d", data))  # '<d' for little-endian double
    return doubles


def plot_histogram(doubles):
    print("we have " + str(len(doubles)) + " data points")
    num_bins = 500
    counts, bin_edges = np.histogram(doubles, bins=num_bins)
    # Compute bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    # normalize the data
    normalized = np.array([float(x) * num_bins / len(doubles) for x in counts])
    # Plot histogram as line graph using matplotlib
    plt.xlabel("Error")
    plt.ylabel("Frequency")
    plt.xlim(0, 40)
    plt.ylim(0, 50)
    plt.plot(bin_centers, normalized)
    plt.show()


# replace 'file.bin' with your binary file
doubles = read_doubles_from_binary_file("debug.data")
plot_histogram(doubles)
print("done!")
