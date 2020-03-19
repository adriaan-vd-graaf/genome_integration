# Calibrating MR-link p values

After a first pass of MR-link and if you have at least 100 and preferably 1,000 uncalibrated p values for different 
genes, it is possible to calibrate them using the script located in `./mr_link/p_value_calibration.py`.

for this you require the `PYMC3` package and the gcc C++ compiler ``g++``. You can install this package using
``` bash
pip3 install pymc3
```
Running a single p value calibration will take up to 30 minutes, but should only have to be performed at the end of 
an analysis, when all the genes are run.


## p value calibration example
After installation of PYMC3 It is possible to run the p value calibration script using the following commands

```shell script
#Run this from the ./mr_link/ directory
python3 p_value_calibration.py --input_file example_files/uncalibrated_p_values_example.txt --output_file calibrated_p_values.txt
```
Which will output calibrated p values in the `calibrated_p_values.txt` file, and accept uncalibrated p values from the
`uncalibrated_p_values_example.txt` file. 

If you want to calibrate _p_ values without computing the beta distribution coefficients, you specify the alpha and beta parameters 
combined with the `--only_calibrate` option in the following way: 
```shell script
#Run this from the ./mr_link/ directory
python3 p_value_calibration.py \
    --input_file example_files/uncalibrated_p_values_example.txt \
    --output_file calibrated_p_values.txt \
    --only_calibrate \
    --alpha_parameter 3.9703 \
    --beta_parameter 0.6106
```
If you use simulated datasets, it's important to calibrate on the null scenarios, then use these parameters on the
non-null scenarios. 




## p_value_calibration.py specifications

The `--input_file` from which p values are calibrated should be a tab separated file with two columns:
1. `exposure_name`
2. `uncalibrated_p_value` 

The `--output_file` is the same file, but with an extra column appended to it:
1. `exposure_name`
2. `uncalibrated_p_value`
3. `calibrated_p_value`
