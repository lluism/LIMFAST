import numpy as np
import os, sys

def get_filepaths(directory, kw):
    """
    This function will generate the file names in a directory
    tree by walking the tree either top-down or bottom-up. For each
    directory in the tree rooted at directory top (including top itself),
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            if kw in filename:
                # Join the two strings in order to form the full filepath.
                filepath = os.path.join(root, filename)
                file_paths.append(filepath)
            else:
                continue

    return file_paths

def write_to_filelist(listpath, listname, kw):

    # Run the above function and store its results in a variable.
    paths = get_filepaths("../Boxes/", kw)
    paths.sort()

    with open(listpath+listname, 'w') as f:
        for i in paths:
            f.write(i)
            f.write('\n')

        f.close()

    return 0

write_to_filelist(listpath='../Redshift_interpolate_filelists/', listname='I_CII158m_256_256Mpc', kw='I_CII158m_')
write_to_filelist(listpath='../Redshift_interpolate_filelists/', listname='I_HI6563A_256_256Mpc', kw='I_HI6563A_')
