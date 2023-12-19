import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from .utils import *
from .petitRADTRANS_init import *

#TODO
# - [ ] Add a function that computes the metallicity and C/O ratio from the initial abundances 


class ACE1D:
    def __init__(self, path_ace):
        """Class that can process the output of an ACE chemistry calculation

        Args:
            path_ace (str): Path to the directory that contains the input and output files of the ACe simulation
        """
        self.path_ace = path_ace
        self.keyword = ""                   
        self.ace_df = None
        self.ace_xr = None
        self.apt_df = None

    def set_path_ace(self, path_ace):
        self.path_ace = path_ace

    def get_path_ace(self):
        return self.path_ace

    def get_ace_df(self):
        return self.ace_df

    def get_ace_xr(self):
        return self.ace_xr

    def load_apt_file(self):
        """Function that loads in the altitude-pressure-temperature file based on the same keyword that is used after loading the .dat-file
        """        
        self.apt_df = load_file_to_dataframe(self.path_ace,
                                            ext = ".apt",
                                            keyword=self.keyword)
    
    def get_apt_file(self):
        return self.apt_df

    def reset_output(self):
        """Function that resets the output variables of this class.
        """
        self.ace_df = None
        self.ace_xr = None
        self.keyword = ""

    def read_1D_ace(self, keyword=""):
        """Loads in the data (.dat) of the ace simulations into the specified formats (pandas dataframe & xarray dataset)

        Args:
            keyword (str, optional): Optional keyword of the .dat-file that is in the ACE director. 
                            This is primarely usefull in the case that there are multiple .dat-files in the ACE directory,
                            such as is the case for different longitudes (e.g. substellar). Defaults to "".

        """

        path = self.get_path_ace()
        self.keyword = keyword

        #Get the paths of all .dat files in the ACE directory
        paths_dat_files = [path + file for file in os.listdir(path) if file.endswith('.dat') and self.keyword in file]
        
        #First, use the first .dat file in the dir to initialize
        #the pressure, altitude and temperature headers as well
        f0 = paths_dat_files[0]
        df = load_one_ACE_dat_file(f0, include_apt_headers=True)

        #Iterate over all other dat files, if there are more than one
        for dat_file in paths_dat_files[1:]:

            #Load in the .dat file without apt headers
            df_temp = load_one_ACE_dat_file(dat_file, include_apt_headers=False)

            #merge togehter with main dataframe
            df = pd.concat([df, df_temp], axis = 1)

        #Change the [] in the variable names to _
        df.rename(columns={'altitude[km]':'altitude_km',
                    'pressure[bar]':'pressure_bar',
                   'Tk[K]':'temperature_K'}, inplace=True)
        
        #store in ace_df
        self.ace_df = df

        #store in ace_xr
        df = df.set_index(['pressure_bar'])
        self.ace_xr = xr.Dataset.from_dataframe(df)


    def plot_TP_profile(self, save = False):
        """Function that plots the PT profile

        Args:
            save (bool, optional): If yes, the figure is saved in the directory of the ACE simulation. Defaults to False.
        """        
        fig, ax = plt.subplots(1, 1, figsize = (10,5))
   
        xr.plot.plot(self.ace_xr["temperature_K"], y = "pressure_bar", yincrease = False, yscale = "log", ax = ax)
        ax.set_xlabel("Temperature")
        ax.set_ylabel("Pressure")

        if save:
            save_path = os.path.join(self.path_ace, "TP_profile_ACE.png")
            fig.savefig(save_path)

    def plot_AP_profile(self, save = False):
        """Function that plots the AP profile

        Args:
            save (bool, optional):  If yes, the figure is saved in the directory of the ACE simulation. Defaults to False.
        """        
        fig, ax = plt.subplots(1, 1, figsize = (10,5))

        xr.plot.plot(self.ace_xr["altitude_km"], y = "pressure_bar", yincrease = False, yscale = "log", ax = ax)

        ax.set_xlabel("Altitude [km]")
        ax.set_ylabel("Pressure [bar]")

        if save:
            save_path = os.path.join(self.path_ace, "AP_profile_ACE.png")
            fig.savefig(save_path)
            
    def get_mean_molecular_weight(self, per_layer = False, molecule_list=None):
        """Function that returns the mean molecular weight (mmw) in amu of a list of molecules.
        WATCH OUT - This is a copy of a function in PAC1D...

        Args:
            per_layer (bool, optional): If true, the mmw will be computed per pressure/altitute layer in the atmosphere. Defaults to False.
            molecule_list (list, optional): List of molecules in string format. If no list is specified, the pac-output is used to determine
             the list of molecules. Defaults to None.

        Returns:
            float: The mean molecular weight
        """
        if molecule_list==None:
            molecule_list = list(self.ace_df.columns)[3:]
        
        # print("molecule_list", molecule_list) 

        total= 0.00
        for mol in molecule_list:

            #Get the weight of this molecule in amu
            mol_weight = get_molecular_weight_molecule(mol)
            
            #get the mean number fraction in the atmosphere
            nf_mol = np.power(10, self.ace_df[mol])

            #Choose the required output...
            if not per_layer:
                nf_mol = nf_mol.mean() #take mean over all dimensions per variable

            #transform to numpy floats
            nf_mol =  np.float64(nf_mol)

            total += mol_weight * nf_mol

        return total
    
    # def get_number_fractions(self):

    #     xrdata = self.xrdat.copy(deep=True)
        
    #     return np.power(10, xrdata.drop(["altitude_km","temperature_K"]))


    # def get_mass_fractions(self):

    #     xr_nf = self.get_number_fractions()

    #     mass_fractions = xr_nf.copy(deep=True)

    #     mu = self.get_mean_molecular_weight(per_layer=False)

    #     for m in self.compute_list_of_molecules():

    #         #Apply formula of pRT to compute mass fractions from number fractions
    #         mu_i = get_molecular_weight_molecule(m, mean=False)

    #         mass_fractions[m].values = (mu_i)/mu * np.float64(xr_nf[m].values)
    # #         print(m,":", mass_fractions[m])
    #     return mass_fractions
            
            
            

#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------

def load_one_ACE_dat_file(path_to_file, include_apt_headers = True):
    """Load in a ACE simulation result (.apt) as a python data frame

    Args:
        path_to_file (str): path to .apt file in the 
        include_apt_headers (bool, optional): If yes, the output dataframe will contain headers. Defaults to True.

    Returns:
        Pandas dataframe: Pandas df that contains the output of the ACE simulation
    """
    f = path_to_file

    #First extract a list of row that start with "!", then select the last row [-1] and cut of the "! " by [2:]
    headers_string = [line for i, line in enumerate(open(f)) if line.startswith('!')][-1][2:]

    # print(headers_string)
    headers_list = headers_string.split()

    #get the indiced of the rows that start with  '!'
    top, end = give_skiprow_top_and_end(f)

    #Extract the data without column names
    max_row = None
    if len(end) != 0:
        max_row = end[0] - len(top)

    #Extract the data without column names
    data = np.loadtxt(f, skiprows= len(top), max_rows=max_row)

    #Include or exclude the altitude, pressure, temperature headers (apt)
    if not include_apt_headers:
        headers_list.remove('altitude[km]')
        headers_list.remove('pressure[bar]')
        headers_list.remove('Tk[K]')
        data = data[:,3:]

    #Make a pandas dataframe
    df = pd.DataFrame(data = data,
                     columns= headers_list,
                     dtype= np.float64)

    return df