import numpy as np
import struct
import os
import gzip
from genome_integration import association
from genome_integration import variants

import genome_integration.association.binary_files

def turn_genetic_association_into_bin(genetic_association):
    try:
        bytes = struct.pack('hiifff',
                            int(genetic_association.chromosome),
                            int(genetic_association.position),
                            float(genetic_association.n_observations),
                            float(genetic_association.beta),
                            float(genetic_association.se),
                            float(genetic_association.r_squared))
        return bytes
    except:
        return b''

def make_binary_file_and_corresponding_bim():

    """
    This makes a binary file which can be read by the binary file class.
    To test how cython would load such a file. randomly putting in files.

    :return: nothing.
    """
    num_snps = 1000

    with gzip.open("tests/test_resources/binary_file_test.dat.gz", "wb") as f:
        f.write(b"testgene\0")
        for i in range(num_snps):
            tmp = association.GeneticAssociation(dependent_name="testgene",
                                                 explanatory_name="1:{}".format(i),
                                                 n_observations=4000,
                                                 beta=np.random.normal(),
                                                 se=np.random.normal(),
                                                 chromosome="1",
                                                 position=i
                                                 )
            f.write(
                turn_genetic_association_into_bin(tmp)
            )

    all_alleles = ['A', 'T', 'C', 'G']
    with open("tests/test_resources/test_bim.bim", "w") as f:
        for i in range(num_snps):
            alleles=np.random.choice(all_alleles, size=2, replace=False)
            f.write("1 1:{} 0 {} {} {}\n".format(i, i, alleles[0], alleles[1]))




def test_reading_binary_file():
    import time

    start = time.time()
    make_binary_file_and_corresponding_bim()
    start = time.time()
    # if this runs the associations are normal.
    bim_file = variants.BimFile("tests/test_resources/test_bim.bim")
    a = genome_integration.association.binary_files.read_bin_file("tests/test_resources/binary_file_test.dat.gz", bim_file)
    end = time.time()
    print(end - start)


if __name__ == '__main__':
    os.chdir("../")
    test_reading_binary_file()