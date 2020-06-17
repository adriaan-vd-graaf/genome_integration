import gzip

"""
This script will turn a 1000g variant file into a .hap and .legend file for hapgen2.
The variant file is an isolated region on chromosome 2 100-105Mb on b37

Only individuals belonging to the european superpopulations excluding finnish individuals are retained
Only bi-alllelic SINGLE nucleotide variants are retained.

To check for sanity the genetic map which hapgen requires is also read, and the positions are checked.
"""

#individuals to keep
non_finnish_europeans = set()
with open("non_finnish_european_individuals.txt", "r") as f:
    for line in f:
        non_finnish_europeans.add(line[:-1])

#the genetic map
genetic_map_positions = []
with open("interpolated_genetic_map.100Mb.to.105Mb.map", "r") as f:
    f.readline()
    for line in f:
        genetic_map_positions.append(line.split()[0])

legend_f = open("europeans_1000g_chr2_100_to_105_Mb.legend",  "w")
legend_f.write("ID pos allele0 allele1\n")

hap_f = open("europeans_1000g_chr2_100_to_105_Mb.hap",  "w")

indices_to_keep = []
i = 0
allele_possibilities = set(["A", "C", "T", "G"])

with gzip.open("ALL.chr2.phase3.100.to.105.Mb.1000g.vcf.gz", "r") as f:
    for line in f:
        if line[0:2] == b"##": #header
            continue
        elif line[0:2] == b"#C": #individual information is here.
            split = line.decode("utf8").split()
            for j in range(len(split)):
                if split[j] in non_finnish_europeans:
                    indices_to_keep.append(j)
        else: #passed all the header lines, now we identify the haplotype information and start printing.
            split = line.decode("utf8").split()

            if split[3] not in allele_possibilities or split[4] not in allele_possibilities:
                i += 1
                continue



            if split[1] != genetic_map_positions[i]: #sanity check.
                raise ValueError(
                    "Positions are not the same between the genetic map and the vcf: gmap {}, vcf {}, line-1 {}".format
                    (genetic_map_positions[i], split[1], i)
                )

            legend_f.write("{} {} {} {}\n".format(split[2], split[1], split[3], split[4]))
            i += 1

            for j in indices_to_keep:
                haps = split[j].split("|")

                # remove copy number variation, will remove at a later point.
                if int(haps[0]) > 1:
                    haps[0] = "1"
                if int(haps[1]) > 1:
                    haps[1] = "1"

                hap_f.write("{} {} ".format(haps[0], haps[1]))

            hap_f.write("\n")