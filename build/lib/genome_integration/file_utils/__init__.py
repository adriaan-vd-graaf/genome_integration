"""
These functions are being used to do some trivial file work..
"""
import gzip

__all__ = ['write_list_to_newline_separated_file',
           'read_newline_separated_file_into_list',
           'read_newline_separated_gz_file_into_list']

__author__      = "Adriaan van der Graaf"
__copyright__   = "Copyright 2017, Adriaan van der Graaf"

def write_list_to_newline_separated_file(in_list, file_name):
    with open(file_name, 'w') as f:
        f.write('\n'.join(in_list))
        f.write('\n') #sometimes R complains, so just adding an extra line.


def read_newline_separated_file_into_list(file_name):
    ret_list = []
    with open(file_name, 'r') as f:
        for line in f:
            ret_list.append(line[:-1])
    return ret_list

def read_newline_separated_gz_file_into_list(filename):
    ret_list = []
    with gzip.open(filename, "rt") as f:
        for line in f:
            ret_list.append(line[:-1])
    return ret_list