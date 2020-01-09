# Calibrating MR-link p values

After a first pass of MR-link and if you have at least 100 and preferably 1,000 uncalibrated p values for different 
genes, it is possible to calibrate them using the script located in `./mr_link/p_value_calibration.py`.

for this you require the `PYMC3` package. You can install this package using
``` bash
pip3 install pymc3
```
Running a single p value calibration will take up to 30 minutes, but should only have to be performed at the end of 
an analysis, when all the genes are run.


## p value calibration example
After installation of PYMC3 It is possible to run the p value calibration script using the following commands

```bash
#Run this from the ./mr_link/ directory
python3 p_value_calibration.py --input_file example_files/uncalibrated_p_values_example.txt --output_file calibrated_p_values.txt
```
Which will output calibrated p values in the `calibrated_p_values.txt` file, and accept uncalibrated p values from the
`uncalibrated_p_values_example.txt` file. 

## p_value_calibration.py specifications

The `--input_file` from which p values are calibrated should be a tab separated file with two columns:
1. `exposure_name`
2. `uncalibrated_p_value` 

The `--output_file` is the same file, but with an extra column appended to it:
1. `exposure_name`
2. `uncalibrated_p_value`
3. `calibrated_p_value`
