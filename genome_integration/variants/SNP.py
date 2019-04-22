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

    def add_minor_allele_frequency(self,  major, minor, freq):
        Adds a minor allele frequency


    _update_alleles(self, snp_data, overwrite)
        updates alleles, but does not update the self.has_allele_data boolean,
        the function update_alleles() will do both.

    def _flip_alleles(self):
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

                if (not snp_data.minor_allele_frequency is None) and \
                        (abs(self.minor_allele_frequency - snp_data.minor_allele_frequency) > 0.5):
                    """
                    This is when the flipped allele frequency is not really right.
                    """
                    warnings.warn("Updated allele frequency difference between snp and reference is really high", RuntimeWarning)

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
