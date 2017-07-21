import argparse
import gzip
import pickle
import os
import math
import time
import shutil

class CohortSize:
    def __init__(self, cohort_loc):
        self.cohort_size = {}
        with open(cohort_loc, "r") as f:
            for line in f:
                split = [x for x in line.split() if x != ""]
                self.cohort_size[split[0]] = int(split[1])

    def find_size(self, cohort_name):
        return self.cohort_size[cohort_name]


def read_in_eqtl_file(frq_loc, eqtl_loc):
    # read in the variants we're interested in.
    variant_list = {}
    with open(frq_loc, "rb") as f:
        f.readline()
        for line in f:
            split = [x for x in line.decode('utf8').split(' ') if x != '']
            variant_list[split[1]] = split

    # now fill the eQTL list.
    eqtl_list = {}
    with gzip.open(eqtl_loc, "rb") as f:
        for line in f:
            split = line.decode('utf8').split('\t')
            try:
                # do this to make sure it's in there.
                a = variant_list[split[1]]
                #now s
                eqtl_list[split[4]] = split
            except:
                pass

    print('Found this many eQTLs in my data:')
    print(str(len(eqtl_list.keys())))
    return (eqtl_list, variant_list)


def find_n_based_on_cohort_line(cohort_line):
    numberdict = {"Cardio": 134,
                 "Bangladeshi": 1404,
                 "CHDWB": 384,
                 "YoungFinns": 1428,
                 "LIFE": 2042,
                 "HVH": 50,
                 "HVN": 50,
                 "Morocco": 175,
                 "SHIP": 955,
                 "SIGN": 115,
                 "Rotterdam": 748,
                 "Fehrmann": 720,
                 "DILGOM": 488,
                 "EGCUT": 818,
                 "EGCUT2": 77,
                 "Sorbs": 513,
                 "BSGS": 329,
                 "KORA": 952,
                 "-": 0
                  }

    #this will throw an exception if it goes wrong, but I expect the cohort to be available always.
    # One other thing, the cohort names which have multiple versions are divided by 2, because I think this may be better,
    # but I could be wrong in thinking that the numbers shown here are of the full cohorts, and the versions are the subset.

    total = 0
    split = cohort_line.split(";")
    for cohort in split:
        total += numberdict[cohort.split("_")[0]]

    return total

def make_probes_(args):

    frq_loc = args.frq
    eqtl_loc = args.eqtl

    #READ overlap with my genotype file. Perhaps nog necessary, but makes your analysis quicker.
    filename_eqtl_list = "eQTL_list.p"
    if not os.path.isfile(filename_eqtl_list):
        eqtl_tuple = read_in_eqtl_file(frq_loc, eqtl_loc)
        pickle.dump(eqtl_tuple, open(filename_eqtl_list,  "wb"))
    else:
        eqtl_tuple = pickle.load(open(filename_eqtl_list, "rb"))

    (eqtl_list, snp_list) = (eqtl_tuple[0], eqtl_tuple[1])

    #MAKE ESD reference file.

    flist_file = open("test_esd.flist", "w")
    flist_file.write("Chr \tProbeID \tGeneticDistance \tProbeBp \tGene \tOrientation \tPathOfEsd\n")

    directory = "probes_from_list"
    probe_list = {}

    for key in eqtl_list.keys():
        eqtl_line = eqtl_list[key]

        probe = eqtl_line[4]

        probe_list[probe] = eqtl_line

        flist_file.write(eqtl_line[2] + " \t")
        flist_file.write("cg" + eqtl_line[4] + " \t")
        flist_file.write("0" + " \t")
        flist_file.write(eqtl_line[6] + " \t")
        flist_file.write("GEN_" + eqtl_line[4] + " \t")
        flist_file.write("-" + " \t")
        flist_file.write(directory + "/" + probe + ".esd\n")

    flist_file.close()

    # Make sure the directory is empty
    shutil.rmtree(directory)
    time.sleep(5)
    os.makedirs(directory)



    # Iterate over all files, and make sure it works.
    writefile_esd = open("placeholder.esd", "w")
    writefile_ma  = open("placeholder.ma", "w")
    with gzip.open(eqtl_loc, "rb") as eqtls:
        for line in eqtls:
            this_eqtl = line.decode('utf8').split('\t')
            if this_eqtl[1] not in snp_list:
                continue
            snpname = this_eqtl[1]
            probe = this_eqtl[4]
            out_filename = directory + "/" + probe\

            if writefile_esd.name != out_filename + ".esd":
                writefile_esd.close()
                writefile_ma.close()
                if os.path.exists(out_filename + ".esd"):
                    writefile_esd = open(out_filename + ".esd", "a")
                    writefile_ma  = open(out_filename + ".ma" , "a")
                else:
                    writefile_esd = open(out_filename + ".esd", "w")
                    writefile_esd.write("Chr\tSNP\tBp\tA1\tA2\tFreq\tBeta\tse\tp\n")

                    writefile_ma  = open(out_filename + ".ma", "w")
                    writefile_ma.write("SNP\tA1\tA2\tfreq\tb\tse\tp\tN\n")

            A2 = snp_list[snpname][3]
            A1 = snp_list[snpname][2]

            snptype_options = [A1 + "/" + A2, A2 + "/" + A1]

            if this_eqtl[8] not in snptype_options:
                print("SNPTYPE did not correspond to SNPfile,found the following:")
                print("\tExpected A1/A2: " + snptype_options[0])
                print("\tFound:          " + this_eqtl[8])
                print("ignoring this snp")
                pass

            freq = float(snp_list[snpname][4])

            z_score = float(this_eqtl[10])

            n = find_n_based_on_cohort_line(this_eqtl[11])

            if A2 != this_eqtl[9]:
                z_score = z_score * -1.0
            # taken from (Zhu et al. 2016)
            beta = z_score / math.sqrt(2 * freq * (1 - freq) * (n + z_score ** 2))
                se = 1 / math.sqrt(2 * freq * (1 - freq) * (n + z_score ** 2))

            writefile_esd.write(this_eqtl[2] + "\t")
            writefile_esd.write(this_eqtl[1] + "\t")
            writefile_esd.write(this_eqtl[3] + "\t")
            writefile_esd.write(A1 + "\t")
            writefile_esd.write(A2 + "\t")
            writefile_esd.write('{:3.5f}'.format(freq) + "\t")
            writefile_esd.write('{:3.5f}\t'.format(beta))
            writefile_esd.write('{:3.5f}\t'.format(se))
            writefile_esd.write(this_eqtl[0] + "\n")

            writefile_ma.write("{}\t{}\t{}\t{:3.5f}\t{:3.5f}\t{:3.5f}\t{}\t{}\n".format(
                this_eqtl[1], A1, A2, freq, beta,se, this_eqtl[0], n)
            )

        writefile_esd.close()
        writefile_ma.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Turn a plink frq file of variants and a FrankeSwertzlab eqtl file into an flist file.')
    parser.add_argument("--frq", required=True, help="plink frq file location")
    parser.add_argument("--eqtl", required=True, help="FrankeSwertzLab eqtl file location")
    args = parser.parse_args()

    main(args)