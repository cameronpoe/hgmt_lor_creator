./hgmt_lor_creator other_images/HGMTDerenzo.phsp b33_effs.csv -e1 -d
python3 plot_histogram.py debug.data Positional\ Error\ \(cm\) First\ Scatter\ Positional\ Error 20 2
./hgmt_debug debug.data 20 40 -hi | tee full_diagnostics/first_scatter_positional_error.txt
./hgmt_lor_creator other_images/HGMTDerenzo.phsp b33_effs.csv -e3 -d
python3 plot_histogram.py debug.data Distance\ \(cm\) First\ Second\ Scatter\ Distance 20 0.2
./hgmt_debug debug.data 20 40 -hi | tee full_diagnostics/first_second_scatter_distance.txt
./hgmt_lor_creator other_images/HGMTDerenzo.phsp b33_effs.csv -e4 -d
python3 plot_histogram.py debug.data Time\ Difference\ \(ns\) First\ Second\ Scatter\ Time\ Difference 0.9 7
./hgmt_lor_creator other_images/HGMTDerenzo.phsp b33_effs.csv -e5 -d
python3 plot_histogram.py debug.data Distance\ \(cm\) First\ Second\ Hit\ Distance 25 0.25
./hgmt_debug debug.data 25 50 -hi | tee full_diagnostics/first_second_hit_distance.txt
./hgmt_lor_creator other_images/HGMTDerenzo.phsp b33_effs.csv -e6 -d
python3 plot_histogram.py debug.data Time\ Difference\ \(ns\) First\ Second\ Hit\ Time\ Difference 6 1.75
./hgmt_debug debug.data 10 40 -hi | tee full_diagnostics/first_second_hit_time_difference.txt
./hgmt_lor_creator other_images/HGMTDerenzo.phsp b33_effs.csv -e7 -d
python3 plot_histogram.py debug.data Energy\ Deposit\ \(KeV\) Energy\ Deposit\ Distribution 511 0.01
./hgmt_debug debug.data 520 52 -hi | tee full_diagnostics/first_scatter_energy_deposit.txt
