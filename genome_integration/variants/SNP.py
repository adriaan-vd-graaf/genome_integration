import warnings

class SNP:
    """
    This is the base SNP class.

    only biallelic variants possible.

    Most of the data  can be missing, as it often is in many datasets that you are trying to parse.


    Attributes
    ----------

    snp_name: str
        Name of the SNP is used for cross comparison if the variant is merged with other snp_info.

    chromosome: int, castable to int or str
        Name of the chromosome, can be a str or an int.

    position: int
        Base pair position on the genome. Internally build b37 is used a lot. So if you're unsure, try and use this.


    major_allele: str
        Used as the more common allele between in the variant (biallelic SNPs only)

    minor_allele: str
        Used as the less common allele between in the variant.

    minor_allle_frequency: float
        float between 0 and 1, inclusive. gives the minor allele frequency.
        By definition, the minor allele should be less often present so should more often than not be <=0.5,
        but alleles can be flipped which would also flip the frequency.

    has_position_data: bool
        if position data is available.

    has_allele_data: bool
        if all alleles are available

    has_frequency_data:
        if frequency data is available.


    Methods
    -------
    add_snp_data(self, snp_data, overwrite=False)
        Adds data from another SNP class-like object.

    add_frequency_data(self, snp_data, flipped, overwrite)
        Adds the frequency data from another SNP object and requires you to say if the alleles are flipped or not.

    add_pos_chr(self, pos, chr):
        adds position information

    update_alleles(self, snp_data, overwrite)
        updates alleles, and updates the self.has_allele_data boolean

    add_minor_allele_frequency(self,  major, minor, freq):
        Adds a minor allele frequency


    _update_alleles(self, snp_data, overwrite)
        updates alleles, but does not update the self.has_allele_data boolean,
        the function update_alleles() will do both.

    _flip_alleles(self):
        Flips alleles -- major becomes minor; minor becomes major. Does not update frequency.

    set_pos_name(self):
        Sets the snp_name attribute to the pos_name which is `{chr}:{position}`


    """

    def __init__(self,
                 snp_name=None,
                 chromosome=None,
                 position=None,
                 major_allele=None,
                 minor_allele=None,
                 minor_allele_frequency=None
                 ):

        self.snp_name = snp_name
        self.chromosome = chromosome
        self.position = position
        self.major_allele = major_allele
        self.minor_allele = minor_allele
        self.minor_allele_frequency = minor_allele_frequency

        self.has_position_data = self.position != None and self.chromosome != None
        self.has_allele_data = self.major_allele != None and self.minor_allele != None
        self.has_frequency_data = self.minor_allele_frequency != None

    def __str__(self):
        return f"{self.snp_name}, {self.chromosome}:{self.position}, major: {self.major_allele}, " \
            f"minor: {self.minor_allele}, MAF: {self.minor_allele_frequency}"

    def _bim_str(self):
        #this only works when using the minor allele as alle
        return f'{self.chromosome}\t{self.snp_name}\t0\t{self.position}\t{self.minor_allele}\t{self.major_allele}'

    def add_snp_data(self, snp_data, overwrite=False):
        """

        This class will return itself with updated snp data.
        It will only change data from a class if the snp_name is the same,
        or if the position is the same.

        Author comment: This is bloody hard to get right.

        :param snp_data, a SNP  object or bigger.:
        :return self:

        """

        # check if we can merge these snps based on names.
        same_name = self.snp_name == snp_data.snp_name
        name_is_position = snp_data.snp_name == "{}:{}".format(self.chromosome, self.position)
        position_is_name = "{}:{}".format(snp_data.chromosome, snp_data.position) == self.snp_name

        if (not same_name) and (not name_is_position) and (not position_is_name):
            raise ValueError("Names of these snps do not match, or the position of these snps do not match")

        # now the snp is accepted, and we look at what we're going to do. in this order.
        position_updated = False


        #add the position data.
        if ((not self.has_position_data) and snp_data.has_position_data) or overwrite:
            self.add_pos_chr(snp_data.position, snp_data.chromosome)
            position_updated = True


        #add the allele data.
        alleles_updated, alleles_flipped = self.update_alleles(snp_data, overwrite)

        #add the frequency data
        frequency_updated = self.add_frequency_data(snp_data, alleles_flipped, overwrite)

        return position_updated, alleles_updated, alleles_flipped, frequency_updated



    def add_frequency_data(self, snp_data, flipped, overwrite):
        """
        updates the frequency data based on a reference.

        :param snp_data:
        :param flipped:
        :return:
        """
        if overwrite:
            self.minor_allele_frequency = snp_data.minor_allele_frequency
            return True

        elif not self.has_frequency_data:
            self.minor_allele_frequency = snp_data.minor_allele_frequency
            self.has_frequency_data = self.minor_allele_frequency != None
            if flipped and self.has_frequency_data:
                self.minor_allele_frequency = 1 - self.minor_allele_frequency

            return True
    
        else:
            if flipped:
                self.minor_allele_frequency = 1 - self.minor_allele_frequency

                if (snp_data.minor_allele_frequency is not None) and \
                        (abs(self.minor_allele_frequency - snp_data.minor_allele_frequency) > 0.5):
                    """
                    This is when the flipped allele frequency is not really right.
                    """
                    warnings.warn(f"High difference in updated allele frequency, "
                                  f"old: {self}, new: {snp_data}",
                                  RuntimeWarning)

            return False


    def add_pos_chr(self, pos, chr):
        """
        adds position information

        :param pos:
        :param chr:
        :return: self
        """

        self.position = pos
        self.chromosome = chr
        self.has_position_data = self.position != None and self.chromosome != None

    def update_alleles(self, snp_data, overwrite):
        """
        updates alleles, and updates the self.has_allele_data boolean
        :param snp_data:
        :param overwrite:
        :return:
        """

        updated, flipped = self._update_alleles(snp_data, overwrite)
        self.has_allele_data = self.major_allele != None and self.minor_allele != None
        return updated, flipped

    def _update_alleles(self, snp_data, overwrite):
        """
        updates alleles, but does not update the self.has_allele_data boolean,
        the function update_alleles() will do both.
        :param snp_data:
        :param overwrite:
        :return:
        """

        updated = False
        flipped = False

        if overwrite:
            self.major_allele = snp_data.major_allele
            self.minor_allele = snp_data.minor_allele
            updated = True
            return updated, flipped

        if (not self.has_allele_data) or self.major_allele == "N" or self.minor_allele == "N" :

            if (self.major_allele is None) and (self.minor_allele is None):

                """
                This is the case where there is no allele data present.
                """

                self.major_allele = snp_data.major_allele
                self.minor_allele = snp_data.minor_allele
                updated = True

                return updated, flipped

            elif not (self.major_allele is None):
                """
                This is the case where only the major allele is present.
                """
                if self.major_allele == snp_data.major_allele:
                    """
                    This is the case where the major allele matches with some reference.
                    """
                    self.minor_allele = snp_data.minor_allele
                    updated = True
                    return updated, flipped
                elif self.major_allele == snp_data.minor_allele:
                    self.minor_allele = snp_data.major_allele
                    self._flip_alleles()
                    updated = True
                    flipped = True
                    return updated, flipped
                else:
                    raise ValueError("Alleles do not match between two snps.")

            elif not (self.minor_allele is None):
                """
                This is the case where only the minor allele is present.
                """
                if self.minor_allele == snp_data.minor_allele:
                    """
                    This is the case where the minor allele matches with the reference allele.
                    """
                    self.major_allele = snp_data.major_allele
                    updated = True
                    return updated, flipped
                elif self.minor_allele == snp_data.major_allele:
                    """
                    This is the case where the minor allele is the opposite of the reference. 
                    """
                    self.major_allele = snp_data.minor_allele
                    self._flip_alleles()
                    updated = True
                    flipped = True
                    return updated, flipped
                else:
                    """
                    This is the case where the alleles do not match with the reference.
                    """
                    raise ValueError("Alleles do not match between two snps.")

        elif self.has_allele_data:
            if (snp_data.major_allele == self.major_allele) and (snp_data.minor_allele == self.minor_allele):
                """
                This is the case where the data in self is the same as the reference.
                """
                return updated,  flipped
            elif (snp_data.major_allele == self.minor_allele) and (snp_data.minor_allele == self.major_allele):
                """
                This is the case where the data in self is the opposite as the reference.
                """
                self._flip_alleles()
                flipped = True
                return updated, flipped
            else:
                """
                This is the case where the data in self does not match the data in the reference.
                """
                raise ValueError("Alleles do not match between two snps.")

        else:
            raise ValueError("This function should never go here. somehow has allele data is None, which is not allowed.")


    def _flip_alleles(self):
        """
        Flips alleles -- major becomes minor, minor becomes major. Does not update frequency.
        :return:
        """
        tmp = self.major_allele
        self.major_allele = self.minor_allele
        self.minor_allele = tmp


    def add_minor_allele_frequency(self,  major, minor, freq):
        """
        Adds a minor allele frequency
        :param major:
        :param minor:
        :param freq:
        :return:
        """

        #if there are no alleles, then just use this.
        if not self.has_allele_data:
            self.minor_allele = minor
            self.major_allele = major
            self.minor_allele_frequency = freq
            self.has_allele_data = True
            self.has_frequency_data = True
            return

        #allele data is present, so we need to check what it is.
        if (self.major_allele == major) and (self.minor_allele == minor):
            self.minor_allele_frequency = freq
            self.has_frequency_data = True
        elif (self.major_allele == minor) and (self.minor_allele == major):
            #need to swap alleles, we just assign them as shown.
            self.major_allele = major
            self.minor_allele = minor

            self.minor_allele_frequency = freq
            self.has_frequency_data = True

        else:
            raise RuntimeError("Alleles do not match in snp" + self.snp_name)

        return

    def set_pos_name(self):
        """
        Sets the snp_name attribute to the pos_name which is `{chr}:{position}`
        :return:
        """
        self.snp_name = "{}:{}".format(self.chromosome, self.position)


class BimFile:
    """
    Reads in a plink formatted bim file that contains many different variants.

    Attributes
    ----------

    bim_results: dict
        Values contains SNP information of all variants in the bim file. keys are the SNP name of the bim file.

    bim_results_by_pos: dict
        Values contains SNP information of all variants in the bim file. keys are the pos name "<chr>:<pos>"
        of the bim file.

    snp_names: list
        ordered list of snp names

    Methods
    -------

    add_frq_information(self, file_name):
        Adds frequency information from a plink file (.frq or .freqx)


    """

    def __init__(self, file_name):

        with open(file_name, 'r') as f:
            all_bim_lines = f.read()

        lines = all_bim_lines.split('\n')
        if lines[-1] == "":
            lines.pop()

        splits_per_line = [[x for x in line.split() if x != ""] for line in lines]


        self.bim_results = {split[1]:
                                        SNP(snp_name=split[1],
                                        chromosome=split[0],
                                        position=split[3],
                                        major_allele=split[5],
                                        minor_allele=split[4],
                                        minor_allele_frequency=None)
                            for split in splits_per_line}
        self.bim_results_by_pos = {f'{split[0]}:{split[3]}':

                                        SNP(snp_name=split[1],
                                        chromosome=split[0],
                                        position=split[3],
                                        major_allele=split[5],
                                        minor_allele=split[4],
                                        minor_allele_frequency=None)
                                   for split in splits_per_line}
        self.snp_names = [split[1] for split in splits_per_line]


    def _write_bim(self, file_name):
        with open(file_name,'w') as f:
            for snp_name in self.snp_names:
                f.write(f'{self.bim_results[snp_name]._bim_str()}\n')


    def add_frq_information(self, file_name):
        """
        Adds frequency information from a plink file (.frq or .freqx)

        :param file_name: file name of the plink frq or frqx file.
        :return: self, with added frq information

        """
        with open(file_name, 'r') as f:
            split = f.readline()[:-1].split("\t")

            # frq file.
            if len(split) == 5:
                for line in f:
                    split = [x for x in line.split() if x != ""]
                    snp_name = split[1]
                    try:
                        self.bim_results[snp_name].add_minor_allele_frequency(split[3], split[2], float(split[4]))
                    except KeyError:
                        try:
                            self.bim_results_by_pos[snp_name].add_minor_allele_frequency(split[3], split[2], float(split[4]))
                        except KeyError:
                            continue

            # frqx file
            elif len(split) == 10:
                for line in f:
                    split = [x for x in line.split() if x != ""]
                    snp_name = split[1]
                    a_one_count = int(split[4])*2 + int(split[5])
                    a_two_count = int(split[6])*2 + int(split[5])
                    if a_one_count <= a_two_count:
                        minor = split[2]
                        major = split[3]
                        maf = float(a_one_count) / float(a_one_count + a_two_count)

                    else:
                        minor = split[3]
                        major = split[2]

                        maf = float(a_two_count) / float((a_one_count + a_two_count))
                    try:
                        self.bim_results[snp_name].add_minor_allele_frequency(major, minor, float(maf))
                    except KeyError:
                        try:
                            self.bim_results_by_pos[snp_name].add_minor_allele_frequency(major, minor, float(maf))
                        except:
                            continue
            else:
                RuntimeError("The frq file header was not in any correct formatting.")