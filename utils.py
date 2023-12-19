import os
import numpy as np
import re
import pandas as pd

def load_file_to_dataframe(path, ext = ".sdf", keyword = ""):
    """This function loads in a .ext file (which is a Pac input file)
    and reads the data into a pandas dataFrame
    To avoid bad usage of the function, I restrict the
    possibilites for the extensions

    Args:
        path (_type_): _description_
        ext (str, optional): _description_. Defaults to ".sdf".
        keyword (str, optional): _description_. Defaults to "".

    Raises:
        Exception: _description_

    Returns:
        _type_: _description_
    """
    
    ext_list = [".apt", ".eddy", ".star"]
    if ext not in ext_list:
        raise Exception(f"E- The requested extension {ext} is not available \n" \
                       f"please select from {ext_list}")

    #search for file extension
    path_file, name_file = get_file_path_and_name(path, ext=ext, keyword=keyword)

    #First extract a list of row that start with "!", then select the last row [-1] and cut of the "! " by [2:]
    headers_string = [line for i, line in enumerate(open(path_file)) if line.startswith('!')][-1][2:]

    # print(headers_string)
    headers_list = headers_string.split()

    #Overwrite with new headers for the .star file
    if ext == ".star":
        headers_list = ['lambda_nm', 'I_nu']

    #get the indiced of the rows that start with  '!'
    top, end = give_skiprow_top_and_end(path_file)

    #Extract the data without column names
    max_row = None
    if len(end) != 0:
        max_row = end[0] - len(top)

    #Extract the data without column names
    data = np.loadtxt(path_file, skiprows= len(top), max_rows=max_row)

    #Make a pandas dataframe
    df = pd.DataFrame(data = data,
                     columns= headers_list,
                     dtype= np.float64)
    return df



def get_file_path_and_name(path, ext=".inp", keyword=""):
    """Function that returns the full path and file name 
    (without .ext) for a given directory input, extension, and keyword.
    If multiple files are found, the first is selected and returned

    Args:
        path (str): path to the parent folder
        ext (str, optional): Extension of end of the file you are looking for. Defaults to ".inp".
        keyword (str, optional): characters that are in the filename you are looking for. Defaults to "".

    Returns:
        file_path: full path of the file you wanted
        name: the file name of the file you wanted
    """
    #search for file extension
    files = [path + file for file in os.listdir(path) if (file.endswith(ext) and keyword in file)]

    #for .inp, we do not want the pac.inp file
    if ext == "inp":
        try:
            files.remove("pac.inp")
        except:
            pass

    if len(files)>1:
        print(f"[WARNING]: The query for files with keyword '{keyword}' and extension '{ext}' has returned {len(files)} results. /nThe first file will be selected")
    file_path = files[0]
    
    name = file_path.split("/")[-1].split(".")[0]
    print(f"[Update]: You have selected {name}{ext}")

    return file_path, name


def get_molecular_weight_molecule(mol_string, compounds = {'He':4.002602,
                                                           'C':12.0110,
                                                           'H':1.0079,
                                                           'O':15.9994,
                                                           'N':14.0067,
                                                           'S':32.065}):
    '''
    Function that returns the (mean) molecular weight in amu. Molecules like "01D" and"03P" are treated seperately. 
    This function can handle molecules starting with a number, such as "2C3H7". It is updated for S-molecules in 2023.

    Args:
        mol_string (str): A string with the name of the molecule e.g. "CH4"

    Returns:
        total (float): The total (mean) molecular weight of the molecule in amu
    '''

    pattern = re.compile(r'([A-Z][a-z]?)(\d+)?')

    #O3P: Atomic oxygen in the ground-level triplet state (highly reactive due to 2 unpaired electrons)
    #O1D: Atomic oxygen in an excited singelt state (rapidly stabilized to O3P by collsion with N2 and O2)
    #Also same comments for N4S and N2D
    if mol_string in ["O3P", "O1D"]:
        return compounds['O']
    if mol_string in ["N4S", "N2D"]:
        return compounds['N']

    total_weight = 0.0
    total_atoms = 0
    for element, n_str in re.findall(pattern, mol_string):

        if element not in ['Y', "Z", 'c']:
            #This list composes some notations that Venot2020 uses and are not physical
            #e.g. for C4H8Y, C2H3CHOZ, cC2H4O
            
            if n_str:
                n = int(n_str)
            else:
                # If there's no number for current element, n_str = '', and the n_str = 1 
                n = 1

            # add to total
            total_weight += compounds[element] * n
            total_atoms += n

    return total_weight

# def count_lines_starting_and_ending_with_exclamation(path):
#     ##USE THIS FUNCTION INSTEAD OF GIVE_SKIPROW_TOP_AND_END - FUTURE WORK
#     start_count = 0
#     end_count = 0
#     with open(path, 'r') as file:
#         lines = file.readlines()
#         for i, line in enumerate(lines):
#             if line.startswith('!'):
#                 if i < len(lines) // 2:
#                     start_count += 1
#                 else:
#                     end_count += 1
#     return start_count, end_count


def give_skiprow_top_and_end(path):
    """
    Function that takes a path to a file and reads that file row per row.
    Every line at the start and at the end of the file, that starts with a "!" is
    counted.

    Parameters
    ----------
    path : str
        path to the file.

    Returns
    -------
    top : list
        list of integers corresponding to the indices of rows starting with "!".
    end : list
        list of integers corresponding to the indices of rows ending with "!".

    """
    #Get all the row indices for the rows that start with !
    line_nb_list = [i for i, line in enumerate(open(path)) if line.startswith('!')]

    if len(line_nb_list) != 1:
        # Split the list in top and end where i is the last index of top
        i = np.argmax(np.diff(line_nb_list))
    else:
        #I always expect there to be a header with '!' so its possible that this is the only line, in that case
        #the argmax function will not work
        i = 0

    top = line_nb_list[:i+1]#; print(top)
    end = line_nb_list[i+1:]#; print(end)

    #Watch out: it is possible that all ! are in the beginning or the end, the gradients are then all 1
    #and so there is no maximum argument
    if len(set(np.diff(line_nb_list)))==1:

        #only !'s in the top of the file
        if line_nb_list[0] ==0:
            top = line_nb_list
            end = []
        #else only '! in the end of the file'
        else:
            top = []
            end = line_nb_list

    return top, end


def ppmol(mol):
    """Function that latex-fies the simple molecule strings, with the intention to print on figures. e.g. by adding subscripts

    Args:
        mol (str): Molecule string

    Returns:
        str: Molecule with subscripts.
    """
    if mol=='NH3':
        return r'NH$_3$'
    if mol=='CH4':
        return r'CH$_4$'
    if mol=='HCN':
        return 'HCN'
    if mol=='C2H2':
        return r'C$_2$H$_2$'
    if mol=='H2O':
        return r'H$_2$O'
    if mol=='H':
        return r'H'
    if mol=='OH':
        return r'OH'
    if mol=='CO2':
        return r'CO$_2$'
    if mol=='CO':
        return r'CO'
    if mol=='NH2':
        return r'NH$_2$'
    if mol=='SO2':
        return r'SO$_2$'
    if mol=='SO':
        return r'SO'
    if mol=='H2S':
        return r'H$_2$S'
    else:
        return mol


def ppside(limb):
    """Function that prints keywords such as "substellar" to a version that can be printed on figures such as "Substellar point"

    Args:
        limb (str): trivial

    Returns:
        str: trivial
    """
    if limb=="substellar":
        return "Substellar point"
    elif limb=="antistellar":
        return "Antistellar point"
    elif limb=="morning":
        return "Morning limb"
    elif limb=="evening":
        return "Evening limb"
    


def read_error_from_datfile(filename):
    """Function that reads the error from the .dat file of an ACE simulation

    Args:
        filename (str): Filename of the .dat file

    Returns:
        float: Error of the ACE simulation file
    """
    with open(filename, 'r') as f:
        # read the contents of the file into a list of lines
        lines = f.readlines()
        # extract the line at index 16 (which is line 17 in 0-indexing)
        line = lines[16]
        # extract the number after the "=" sign
        number = line.split("=")[1].strip()
        # convert the number to a float and return it
        return float(number)
    
def give_errors(list_files):
    """Function that reads the errors from a list of .dat files

    Args:
        list_files (list): List of paths to .dat files

    Returns:
        list: List of errors corresponding to the entries in list_files
    """
    res = []
    for file in list_files:
        res.append(read_error_from_datfile(file))
    return res

def find_datfiles_and_sort(directory):
    """Function that finds all the .dat files in a directory and sorts them based on the integer value in the filename

    Args:
        directory (str): Directory in which the .dat files are stored

    Returns:
        list: List of paths to the .dat files
    """

    files = []
    for filename in os.listdir(directory):
        if filename.endswith("_1.dat"):
            files.append(os.path.join(directory, filename))
    # sort the list of files based on the integer value in the filename
    files.sort(key=lambda x: int(x.split("_")[-2]))
    return files


def roll_to_degrees(xr_array):
    """This function transforms the longitude coordinates from 0 to 180 to  -90 to 90 (or reverse?)

    Args:
        xr_array (_type_): _description_
    """
    
    half = int(xr_array.sizes["longitude_npi"]/2)

    #First roll the variables, without altering the coordinates
    tmp = xr_array.roll(longitude_npi = half, roll_coords=False)

    #Remap the dimension coordinates
    remapped_lon =  (((xr_array.longitude_npi*180 + 180) % 360) - 180)
    rolled_long = remapped_lon.roll(longitude_npi = half, roll_coords=True)
    tmp = tmp.assign_coords(longitude_npi = rolled_long)

    return tmp

def get_subset_around_long(xr_array, center = 0, width=0.1):
    """Function that slices an xarray dataset around a certain longitude

    Args:
        xr_array (xarray.DataArray): Data structure of the xarray dataset
        center (int, optional): Central longitude. Defaults to 0.
        width (float, optional): Width of the slice around the central longitude. Defaults to 0.1.

    Returns:
        xarray.DataArray: DataArray that is sliced around the central longitude
    """
    condition = False

    ref = center #This is the longitude [0,2] which you want to look at
    if width:
        hw = width/2
        
        if 0<(ref-hw)<2 and 0<(ref+hw)<2:
            condition = (xr_array.longitude_npi > ref-hw) & (xr_array.longitude_npi < ref+hw)
        elif (ref-hw)<0 and 0<(ref+hw)<2:
            condition = (xr_array.longitude_npi < ref+hw) | (xr_array.longitude_npi > 2-ref-hw)
        elif 2>(ref-hw)>0 and (ref+hw)>2:
            condition = (xr_array.longitude_npi > ref-hw) & (xr_array.longitude_npi < (ref+hw)-2)
    
    #Check if condition is not only Falses
    if (not np.any(condition)) or (width is None):
        xr_array_sliced = xr_array.sel(longitude_npi = center , method="nearest")
    else:
        xr_array_sliced = xr_array.sel(longitude_npi=xr_array.longitude_npi[condition])
        xr_array_sliced=xr_array_sliced.mean(dim = "longitude_npi")
    
    return xr_array_sliced