import argparse
import numpy as np
import pymc3 as pm
import scipy.stats

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_file',
                        type=str,
                        help="File name with uncalibrated p values",
                        required=True)

    parser.add_argument("--output_file",
                        type=str,
                        help="the file name where the MR-link results will be output.",
                        required=True)

    parser.add_argument("--tuning_steps",
                        type=int,
                        help="The number of tuning steps in the MCMC inference"
                             "Increase this if you get an acceptance probability warning.",
                        default=1000000)

    parser.add_argument("--only_calibrate",
                        action='store_true',
                        help="Use this, if you already have alpha and beta values computed")

    parser.add_argument("--alpha_parameter",
                        type=float,
                        default=None,
                        help="alpha parameter of the beta distribution, from which to calibrate p values. "
                             "Only possible if the parameter has been precomputed and --only_calibrate option "
                             "is specified.")

    parser.add_argument("--beta_parameter",
                        type=float,
                        default=None,
                        help="beta parameter of the beta distribution, from which to calibrate p values. "
                             "Only possible if the parameter has been precomputed and --only_calibrate option "
                             "is specified.")

    args = parser.parse_args()

    ##runtime checks
    if (args.alpha_parameter is not None or args.beta_parameter is not None) and not args.only_calibrate:
        raise ValueError("Alpha and beta parameters should not be specified if --only_calibrate is not specified.")

    if args.only_calibrate and not isinstance(args.alpha_parameter, float):
        raise ValueError("If --only_calibrate is specified, --alpha_parameter needs to be specified too")

    if args.only_calibrate and not isinstance(args.beta_parameter, float):
        raise ValueError("If --only_calibrate is specified, --beta_parameter needs to be specified too")



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

    if args.only_calibrate:
        posterior_alpha = args.alpha_parameter
        posterior_beta = args.beta_parameter
        print(f"Using pre-specified alpha and beta parameters: alpha: {posterior_alpha} and beta: {posterior_beta}")
    else:
        model = pm.Model()

        with model:
            alpha = pm.Uniform("alpha", lower=0, upper=10)
            beta = pm.Uniform("beta", lower=0, upper=10)

            y_like = pm.Beta('distribution', alpha=alpha, beta=beta, observed=uncalibrated_p_values)

            # draw 100,000 posterior samples after 1,000,000 tuning steps.

            trace = pm.sample(100000, tune=1000000, chains=2)

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

    print(f"Finished analysis on {args.input_file}, with a posterior alpha: {posterior_alpha} and beta: {posterior_beta}")