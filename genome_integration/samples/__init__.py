__author__      = "Adriaan van der Graaf"


class Sample:
    """
    Implements a sample. only name is stored.

    Attributes
    ----------

    Name: str sample name

    Phenotype: can be a value, an array or dictionary of values, but initialized as None

    """

    def __init__(self, name, phenotype=None):

        self.name = name
        self.phenotype = phenotype

