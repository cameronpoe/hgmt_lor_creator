echo "Running Diagnostics"
./hgmt_lor_creator data/HGMTDerenzo.phsp kapton_effs.csv -e0 -e1 -e2 -d | tee full_diagnostics/detector_diagnostics.txt
./hgmt_debug data/debug0.data 13 13 -hi | tee full_diagnostics/scatter_detector.txt
python3 plot_bars.py full_diagnostics/scatter_detector.txt Detector\ ID Scatters\ Occured
./hgmt_debug data/debug1.data 13 13 -hi | tee full_diagnostics/first_scatter_detector.txt
python3 plot_bars.py full_diagnostics/first_scatter_detector.txt Detector\ ID First\ Scatters\ Occured
python3 plot_histogram.py data/debug2.data Lor\ Center\ Error\ \(cm\) First\ Scatter\ Positional\ Error 10 2
./hgmt_debug data/debug2.data 20 1000 -hi | tee full_diagnostics/lor_center_error.txt
echo "All Tasks Complete"
