echo "Running Diagnostics"
./hgmt_lor_creator data/HGMTDerenzo.phsp kapton_effs.csv -e0 -e1 -d | tee full_diagnostics/detector_diagnostics.txt
./hgmt_debug data/debug0.data 13 13 -hi | tee full_diagnostics/scatter_detector.txt
python3 plot_bars.py full_diagnostics/scatter_detector.txt Detector\ ID Scatters\ Occured
./hgmt_debug data/debug1.data 13 13 -hi | tee full_diagnostics/first_scatter_detector.txt
python3 plot_bars.py full_diagnostics/first_scatter_detector.txt Detector\ ID First\ Scatters\ Occured
echo "All Tasks Complete"
