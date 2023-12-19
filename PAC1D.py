import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import numpy as np
from .utils import *
from .petitRADTRANS_init import *
from .PAC_super import PAC_super
import astropy.constants as cte

class PAC1D(PAC_super):

    def __init__(self, outputdir, path_figsave=None):
        """Class that processes the output of a 1D PAC chemistry simulation

        Args:
            outputdir (str): Path to the directory containing all the input and outputfiles of the PAC simulation
            path_figsave (str): Path to a directory where the generated figures can be saved.
        """    

        #Init the super class
        super().__init__(outputdir, path_figsave)

        self.converged_output_nb, self.outputfile_path = self._find_converged_output_file(outputdir)
        self.resolution = None


    def get_pac1d_df(self):
        return self.pac1d_df

    def get_converged_output_nb(self):
        return self.converged_output_nb

    def get_outputfile_path(self):
        return self.outputfile_path
    
    def reset_outputfile(self):
        """Function that resets the outputfile to the converged outputfile
        """
        self.converged_output_nb, self.outputfile_path = self._find_converged_output_file(self.outputdir)
    
    def set_converged_output_nb(self, number=None):
        """This function allows you to overwrite the converged_output_nb to a string of your choice

        Args:
            number (_type_, optional): _description_. Defaults to None.
        """
        if number==None:
            self.converged_output_nb, self.outputfile_path = self._find_converged_output_file(self.outputdir)
            print(f"WARNING: No number is specified bu the user, and so we choose the converge outputfile with number {self.converged_output_nb}")

        else:
            self.converged_output_nb = number

        # Don't forget to set the path as well
        self.outputfile_path, _ = get_file_path_and_name(path=self.outputdir,
                                            ext=self.converged_output_nb + "_1.dat",
                                            keyword="out"
                                            )

    def read_1D_pac(self):
        """
        Read a 1-dimensional PAC simulation output file, and sets two class variables.
        Function adopted from Robin Baeyens.
        """

        # Parse header
        my_attrs = {}
        with open(self.outputfile_path, 'r') as f:
            for line in f:
                items = line.split()
                if items[0] == '!':
                    if len(items) > 1:
                        if items[1] == 'species': my_attrs['specfile'] = items[4]
                        elif items[1] == '(z,p,T)': my_attrs['zabfile'] = items[8]
                        elif items[1] == 'reaction': my_attrs['reacfile'] = items[8]
                        elif items[1] == 'stellar': my_attrs['starfile'] = items[5]
                        elif items[1] == 'photo': my_attrs['photodiss-list'] = items[8]
                        elif items[1] == 'Star-planet': my_attrs['d (AU)'] = items[4]
                        elif items[1] == 'Star': my_attrs['R_star (Rsun)'] = items[4]
                        elif items[1] == 'Planet':
                            if items[2] == 'radius': my_attrs['R_p (Rearth)'] = items[4]
                            elif items[2] == 'mass': my_attrs['M_p (Mearth)'] = items[4]
                        elif items[1] == 'Zenith': my_attrs['zenith angle'] = items[4]
                        elif items[1] == 'Diffusion': my_attrs['vertical mixing'] = items[5]
                        elif items[1] == 'Photochemistry': my_attrs['photochemistry'] = items[5]
                        elif items[1] == 'Numerical': my_attrs['numerical method'] = items[5]
                        elif items[1] == 'Continuity': my_attrs['integration time'] = items[7]
                        elif items[1] == 'with': my_attrs['relative error'] = items[7]
                        elif items[1] == 'Steady': my_attrs['automatic steady state'] = items[6]
                        elif items[1] == 'height[km]': header_names = items[1:]
                else:
                    break

        # Read the remainder of the file with the header names known
        df = pd.read_csv(self.outputfile_path, comment='!', delim_whitespace=True,
                        names=header_names)
                        
        df = df.sort_values(by="pressure[bar]")

        #Change the [] in the variable names to _
        df.rename(columns={'height[km]':'height_km',
                    'pressure[bar]':'pressure_bar',
                   'temperature[K]':'temperature_K'}, inplace=True)

        self.pac1d_df = df
        self.my_attrs = my_attrs
        self.R_pl = float(self.my_attrs["R_p (Rearth)"]) * cte.R_earth.cgs.value
        self.M_pl = float(self.my_attrs["M_p (Mearth)"]) * cte.M_earth.cgs.value
        self.R_star = float(self.my_attrs["R_star (Rsun)"]) * cte.R_sun.cgs.value
        self.gravity = self._compute_gravity()

        self.temperatures = self.pac1d_df["temperature_K"].to_numpy()
        self.pressures = self.pac1d_df["pressure_bar"].to_numpy()
        self.resolution = len(self.pressures)


    def plot_convergence(self, molecules=None):
        """Function that plots the convergence evolution of the PAC simulation 
        """
        #prepare the list of files
        list_files = find_datfiles_and_sort(self.outputdir)
        cmap = plt.get_cmap("viridis")
        colors = list(cmap(np.linspace(0.01,1, len(list_files)+1)))       
        c_iter = iter(colors)

        #Init plot
        if molecules!=None:
            fig, [ax, ax2] = plt.subplots(1,2)
        else:
            fig, ax = plt.subplots(1,1)

        # Plot the relative error
        list_errors = give_errors(list_files)
        xs = np.linspace(1, len(list_errors), num=len(list_errors))
        for x, err in zip(xs, list_errors):
            ax.semilogy(x, err, ls="", ms=10, marker="o",
                    color=next(c_iter))
            
        ax.set_ylim([0.5*min(list_errors), 1e5])
        ax.set_xlabel("file number")
        ax.set_ylabel("relative error")
        ax.set_yscale('log')
        ax.grid(True)


        #Plot the molecules
        if molecules!=None: 
            c_iter = iter(colors)

            label_flag = True
            for file in list_files:
                c = next(c_iter)
                self.outputfile_path = file
                self.read_1D_pac()
                df = self.pac1d_df
                
                #only add a label to the first one
                for mol in molecules:
                    
                    label=None
                    if label_flag==True:
                        label=mol
                    try:
        #                 print(file)
                        ax2.semilogy(df[mol], df["pressure_bar"],
                                    c = c, label = label, 
                                    lw=2.5)
                    except:
                        print(f"Something went wrong with {mol}")
                        pass
                    
                label_flag=False
                
            ax2.set_ylabel(r"")
            ax2.set_xlabel(r"")        
            ax2.set_ylim([200, 1e-8])
            ax2.set_title('')
            ax2.set_xlabel("Log molar fraction")
            ax2.set_ylabel(r"Pressure [bar]")

            #Reset the outputfile
            self.reset_outputfile()

    def _find_converged_output_file(self, outputdir):
        """This function computes the file number at which the pac simulation is converged by looking at the minimal value for the relative
        error in the second half of outputs
        
        Args:
            outputdir (str): Path to the directory that contains all the outputfiles

        Returns:
            str: Two-digit number in string format that represents the number of the highest outputfile
            str: Path to the highest outputfile
        """
        
        list_files = find_datfiles_and_sort(outputdir)
        list_errors = give_errors(list_files)

        #only keep the second half of files.
        N = len(list_errors)
        list_files = list_files[int(N/2):]
        list_errors = list_errors[int(N/2):]
        i_lowest = np.argmin(np.array(list_errors))

        #Select
        outputfile_path = list_files[i_lowest]
        converged_output_nb = outputfile_path.split("_")[-2]
            
        return converged_output_nb, outputfile_path


    def _find_highest_numbered_file(self, path):
        """Function that finds the file with the last outputfile in the directory.
        Only to be used in case that the model has converged to the lowest relative error at the end.
        This rarely the case, so use another function called _find_converged_output_file()

        Args:
            path (str): Path the parentfolder, containing the .dat output files of the PAC code (outputdir)

        Returns:
            str: Two-digit number in string format that represents the number of the highest outputfile
            str: Path to the highest outputfile
        """

        pattern = os.path.join(path, "*_1.dat")
        files = glob.glob(pattern)
        if not files:
            return None
        highest_number = max(int(os.path.basename(file).split('_')[-2]) for file in files)
        highest_file = next((file for file in files if int(os.path.basename(file).split('_')[-2]) == highest_number), None)

        return f"{highest_number:02}" , highest_file


    def get_mean_molecular_weight(self, per_layer = False, molecule_list=None):
        """Function that returns the mean molecular weight (mmw) in amu of a list of molecules.

        Args:
            per_layer (bool, optional): If true, the mmw will be computed per pressure/altitute layer in the atmosphere. Defaults to False.
            molecule_list (list, optional): List of molecules in string format. If no list is specified, the pac-output is used to determine
             the list of molecules. Defaults to None.

        Returns:
            float: The mean molecular weight
        """

        if molecule_list==None:
            molecule_list = list(self.pac1d_df.columns)[3:]

        total= 0.00
        for mol in molecule_list:

            #Get the weight of this molecule in amu
            mol_weight = get_molecular_weight_molecule(mol)
            
            #get the mean number fraction in the atmosphere
            nf_mol = np.power(10, self.pac1d_df[mol])

            #Choose the required output...
            if not per_layer:
                nf_mol = nf_mol.mean() #take mean over all dimensions per variable

            #transform to numpy floats
            nf_mol =  np.float64(nf_mol)

            total += mol_weight * nf_mol

        return total



    def plot_chemistry(self,
                list_molecules = ['C2H2', 'CO2', 'OH', 'HCN', 'NH3', 'CH4', 'CO', 'H2O', 'H2', "SO2", "H2S"],
                ACE = None,
                save=False):
        """Function that plots the basic chemistry output of the PAC simulation.

        Args:
            list_molecules (list, optional): List of molecules (string) that should be plotted on the figure.
                                             Defaults to ['C2H2', 'CO2', 'OH', 'HCN', 'NH3', 'CH4', 'CO', 'H2O', 'H2', "SO2"].
            save (bool, optional): If true, the figure is save in the directory specified by the init. Defaults to False.
        """
        
        #Init plot
        fig, axs = plt.subplots(1,2,figsize=(12,6),
                                sharey=True, sharex=False,
                                gridspec_kw = {'wspace':0.0, "width_ratios":[1, 2]}
                            )       
        axs = axs.reshape(-1)
        
        #Init colors
        cmap = plt.get_cmap("tab20")
        colors = list(cmap(np.linspace(0.01,1, len(list_molecules)+1)))       
        c_iter = iter(colors)
        
        
        #Plot the T-P profile
        ax = axs[0]
        df_pac1d = self.get_pac1d_df()
        ps = df_pac1d["pressure_bar"]
        Ts = df_pac1d["temperature_K"]
        ax.semilogy(Ts, ps)
        ax.grid(alpha=0.5)
        
        #Plot the molecules
        ax = axs[1]
        if ACE!=None:
            ace_df = ACE.get_ace_df()

        for m in list_molecules:

            #Try plotting
            try:
                c = next(c_iter)
                ax.semilogy(df_pac1d[m], df_pac1d["pressure_bar"],
                                c = c, label = m, lw=2.5)
            except:
                print(f"Something went wrong with {m} for PAC")

            if ACE!=None:
                try:
                    ax.semilogy(ace_df[m], ace_df["pressure_bar"],
                                    ls="--",c = c, lw=1.5)
                except:
                    print(f"Something went wrong with {m} for ACE")

        c_iter = iter(colors)

        #Fix the legends
        legend = axs[1].legend(loc="lower right", labelcolor='linecolor',
                        handletextpad=0, handlelength=0, 
                        frameon=False)
        for item in legend.legendHandles:
            item.set_visible(False)
                        
                
        #Deal with the axis
        axs[1].set_xlabel("Log molar fraction")
        axs[0].set_xlabel("Temperature [K]")
        axs[0].set_ylabel(r"Pressure [bar]")
        axs[1].set_xlim([-15, 0])
        axs[0].set_ylim([np.max(df_pac1d["pressure_bar"]), np.min(df_pac1d["pressure_bar"])])
        # axs[0].set_xlim([300, 3500])
        # axs[0].set_xticks([1000, 2000, 3000])

        if save and self.path_figsave!=None:
            fig.savefig(self.path_figsave + "Chemistry_plot" + ".png",
                    bbox_inches='tight')


#            ____  ______
#     ____  / __ \/_  __/
#    / __ \/ /_/ / / /   
#   / /_/ / _, _/ / /    
#  / .___/_/ |_| /_/     
# /_/                    for petitRADTRANS usage and mass fractions


    def set_mass_fractions(self, PRT_ATM, additional_massFractions=None):
        """Function that computes the mass fractions of the molecules in the PAC simulation, 
        and stores them in the class variable mass_fractions

        Args:
            PRT_ATM (apu.PRT_ATM object): Object that controls the petitRADTRANS atmosphere
            additional_massFractions (dict, optional): Dictionary containing additional mass fractions that should be included. Defaults to None.

        Raises:
            Exception: Raises an exception if the mass fractions of some molecules are not included in the PAC simulation and not specified in the additional_massFractions dictionary.
        """

        #Mean molecular weight
        mmw = self.get_mean_molecular_weight(per_layer=True)

        mass_fractions = {}
        mass_fractions_missing = []

        #Species for which mass fractions should be computed
        mf_molecules = PRT_ATM.give_mf_molecules()
        
        for m in mf_molecules:

            try:
                
                #Apply formula of pRT to compute mass fractions from number fractions
                mu_i = get_molecular_weight_molecule(m)

                #Read the number fraction in log-format
                nf = self.pac1d_df[m].to_numpy()
                
                #Spit out the correct pRT name for the molecule e.g. SO2 (pac) becomes SO2_ExoAmes (pRT)
                pRT_key = translate_PAC_to_pRT_molecule(m)
                
                #Add to mass fraction dictionary
                mass_fractions[pRT_key] = (mu_i)/mmw * np.float64(np.power(10, nf))

            except:
                mass_fractions_missing.append(m)
            
        # print("The mass fractions that should be additionally included are:", mass_fractions_missing)

        if additional_massFractions!=None:
            for m in list(additional_massFractions.keys()):            
                mass_fractions[m] = additional_massFractions[m]
                mass_fractions_missing.remove(m)

        if len(mass_fractions_missing)!=0:
            raise Exception("The following molecules are not included in the mass fractions:", mass_fractions_missing)

        self.mass_fractions = mass_fractions
        self.PRT_ATM = PRT_ATM
        self.MMW = mmw