import os
import sys
import subprocess


def test_mr_link_results_no_causal_effect():


    test_dir = sys.path[0]
    mr_link_path = test_dir + "/../mr_link/"

    output_file  = f"{test_dir}/temp_data/no_causal_effect_example.txt"
    subprocess.run( #this example should run correctly
        ['python3', f'{mr_link_path}/MRlink.py',
         "--outcome_bed_file", f"{mr_link_path}example_genotypes/outcome_cohort",
         "--reference_bed", f"{mr_link_path}example_genotypes/reference_cohort",
         "--exposure_summary_statistics", f"{mr_link_path}example_files/no_causal_effect_exposure_sumstats.txt",
         "--outcome_phenotype_file", f"{mr_link_path}example_files/no_causal_effect_outcome_pheno.txt",
         "--temporary_location_prepend", f"{test_dir}/temp_data/tmp_loc",
         "--p_val_iv_selection_threshold", "5e-8",
         "--output_file", output_file,
         "--ensg_id", "ENSG00000000000"], check=True
    )

    correct_header = "ensembl_name\tmethod\tbeta\tse\tp_value\tn_ivs\tiv_summary\n"
    reference_result = 'ENSG00000000000\tMR-link_uncalibrated\t-0.01278\t0.10784\t9.057e-01\t3\trs58572539;C;-0.40832;0.02737;0.163,rs77150110;A;0.48895;0.05324;0.037,rs3771166;A;-0.26576;0.02074;0.393\n'

    with open(output_file, 'r') as f:
        header = f.readline()
        assert header == correct_header
        result = f.readline()
        assert result == reference_result


    os.remove(output_file)


def test_mr_link_results_yes_causal_effect():
    test_dir = sys.path[0]
    mr_link_path = test_dir + "/../mr_link/"

    output_file = f"{test_dir}/temp_data/yes_causal_effect_example.txt"
    print(output_file)
    subprocess.run(  # this example should run correctly
        ['python3', f'{mr_link_path}/MRlink.py',
         "--outcome_bed_file", f"{mr_link_path}example_genotypes/outcome_cohort",
         "--reference_bed", f"{mr_link_path}example_genotypes/reference_cohort",
         "--exposure_summary_statistics", f"{mr_link_path}example_files/yes_causal_effect_exposure_sumstats.txt",
         "--outcome_phenotype_file", f"{mr_link_path}example_files/yes_causal_effect_outcome_pheno.txt",
         "--temporary_location_prepend", f"{test_dir}/temp_data/tmp_loc",
         "--p_val_iv_selection_threshold", "5e-8",
         "--output_file", output_file,
         "--ensg_id", "ENSG00000000000"], check=True
    )

    correct_header = "ensembl_name\tmethod\tbeta\tse\tp_value\tn_ivs\tiv_summary\n"
    reference_result = 'ENSG00000000000\tMR-link_uncalibrated\t0.41500\t0.14734\t4.852e-03\t2\trs11679856;G;0.31542;0.02761;0.197,rs62156384;G;-0.73949;0.02757;0.197\n'

    with open(output_file, 'r') as f:
        header = f.readline()
        assert header == correct_header
        result = f.readline()
        assert result == reference_result

    os.remove(output_file)


test_mr_link_results_yes_causal_effect()
test_mr_link_results_no_causal_effect()
