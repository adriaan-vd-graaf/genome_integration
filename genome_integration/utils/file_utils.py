import gzip

def write_list_to_newline_separated_file(in_list, file_name):
    """
    Convenience function: 
    Writes list into newline separated file.
    
    :param in_list: list of str
        list to write away. Lines need to be not newline separated.
    :param file_name: str 
        out file name 
    :return: None
    """
    with open(file_name, 'w') as f:
        f.write('\n'.join(in_list))
        f.write('\n') # sometimes R complains, so just adding an extra line, for good measure.


def read_newline_separated_file_into_list(file_name):
    """
    Reads a file into a list of strings.
    :param file_name: str
        filename to read    
    :return: list of str
        newline characters are removed.
    """
    ret_list = []
    with open(file_name, 'r') as f:
        for line in f:
            ret_list.append(line[:-1])
    return ret_list

def read_newline_separated_gz_file_into_list(filename):
    """
    Reads a gzipped file into a list of strings.
    
    :param file_name: str
        filename to read    
    :return: list of str
        newline characters are removed.
    """
    
    ret_list = []
    with gzip.open(filename, "rt") as f:
        for line in f:
            ret_list.append(line[:-1])
    return ret_list