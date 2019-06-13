import argparse
import numpy as np
import pymc3 as pm
import scipy.stats

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_file',
                        type=str,
                        help="File name with uncalibrated p values")

    parser.add_argument("--output_file",
                        type=str,
                        help="the file name where the MR-link results will be output."
                        )

    args = parser.parse_args()

    gene_names = []
    uncalibrated_p_values = []

    correct_input_header = 'exposure_name\tuncalibrated_p_value\n'
    with open(args.input_file, 'r') as f:
        header = f.readline()

        if header != correct_input_header:
            raise ValueError(
                "The header of --input_file should be exactly the following: \'exposure_name\\tuncalibrated_p_value\\n\'")

        for line in f:
            split = line.split()
            gene_names.append(split[0])
            uncalibrated_p_values.append(float(split[1]))
            if not (1.0 > float(split[1]) > 0.0):
                raise ValueError("uncalibrated P values should be between 0 and 1 inclusive.")

    uncalibrated_p_values = np.asarray(uncalibrated_p_values, dtype=float)

    # values of exactly 1.0 are
    uncalibrated_p_values[uncalibrated_p_values == 1.0] = np.nextafter(1.0, -2) #


    model = pm.Model()

    with model:
        alpha = pm.Uniform("alpha", lower=0, upper=10)
        beta = pm.Uniform("beta", lower=0, upper=10)

        y_like = pm.Beta('distribution', alpha=alpha, beta=beta, observed=uncalibrated_p_values)

        # draw 100,000 posterior samples after 1,000,000 tuning steps.
        step = pm.Slice()
        trace = pm.sample(100000, slice=step, tune=1000000, chains=2)

    result = pm.summary(trace)

    posterior_alpha = result['mean'][0]
    posterior_beta = result['mean'][1]

    print(f"Finished fitting beta distribution, posteriors: alpha {posterior_alpha:.4f}, beta {posterior_beta:.4f}")
    calibrated_p_values = scipy.stats.beta.cdf(uncalibrated_p_values, a=posterior_alpha, b=posterior_beta)
    out_header = correct_input_header = 'exposure_name\tuncalibrated_p_value\tcalibrated_p_value\n'
    with open(args.output_file, 'w') as f:
        f.write(out_header)
        for name, uncalibrated, calibrated in zip(gene_names, uncalibrated_p_values, calibrated_p_values):
            f.write(f"{name}\t{uncalibrated}\t{calibrated}\n")