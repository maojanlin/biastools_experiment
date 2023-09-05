
#bias_report="./bias_report/wgs.vg_pop5.chr20.bias"
#bias_report="./bias_report_octopus/wgs.octopus.vg_major.chr20.bias"
#bias_report="./bias_report/wgs.giraffe_1KGP001.chr20.bias"
#bias_report="./bias_report_octopus/wgs.octopus.giraffe_1KGP001.chr20.bias"
#bias_report="./bias_report/wgs.bt2_local.chr20.bias"
#bias_report="./bias_report_adaptive/wgs.giraffe_1KGP001_all.chr20.adaptive.bias"
#bias_report="./bias_report_overlap/wgs.giraffe_1KGP001.chr20.overlap.bias"
#bias_report="./bias_report_bothinvariant/wgs.bt2.chr20.overlap.round.bias"
bias_report="./bias_report_bothinvariant/wgs.bt2.chr20.bothinvariant.bias"
#bias_report="./bias_report_overlap/wgs.bwa_e2e.chr20.bias"
#bias_report="./bias_report_naive/wgs.bt2.chr20.naive.bias"

echo "[BIASTOOLS] Plotting"
# making the report of different categories of bias
python3 ../biastools/golden_graph_report.2.py   -mb ${bias_report}".SNP"
python3 ../biastools/golden_graph_report.2.py   -mb ${bias_report}".gap"
# plot the measures with NMB and NAB
python3 ../biastools/golden_graph.py          -mb ${bias_report}".SNP"
python3 ../biastools/golden_graph.py          -mb ${bias_report}".gap"
# generate the balance files without golden / show the mapQ and balance difference between bias sites if -SIM is specified
python3 ../biastools/depth_balance_grapher.py -mb ${bias_report}".SNP" -sim
python3 ../biastools/depth_balance_grapher.py -mb ${bias_report}".gap" -sim

