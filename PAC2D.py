import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
from .utils import *
from .petitRADTRANS_init import *
from .PAC_super import PAC_super
import astropy.constants as cte


class PAC2D(PAC_super):
    def __init__(self, outputdir, path_figsave=None):

        # Init the super class
        super().__init__(outputdir, path_figsave)
        self.pac2d_xr = None

        #Other variables
        self.lpt_xr = None
        self.mws = None #mean wind speed

        #Resolution of the PAC2D simulation
        self.resolution_z = None
        self.resolution_l = None
        self.rotation_period = None

    def get_resolution_z(self):
        return self.resolution_z
    def get_resolution_l(self):
        return self.resolution_l
    
    def _read_my_attrs(self, f0):

        # Parse header
        my_attrs = {}
        with open(f0, 'r') as f:
            for line in f:
                items = line.split()
                if items[0] == '!':
                    if len(items) > 1:
                        if items[1] == 'species': my_attrs['specfile'] = items[4]
                        elif items[1] == '(z,p,T)': my_attrs['zabfile'] = items[8]
                        elif items[1] == '(l,p,T)': my_attrs['lptfile'] = items[6]
                        elif items[1] == 'reaction': my_attrs['reacfile'] = items[8]
                        elif items[1] == 'stellar': my_attrs['starfile'] = items[5]
                        elif items[1] == 'photo': my_attrs['photodiss-list'] = items[8]
                        elif items[1] == 'Star-planet': my_attrs['d (AU)'] = items[4]
                        elif items[1] == 'Star': my_attrs['R_star (Rsun)'] = items[4]
                        elif items[1] == 'Planet':
                            if items[2] == 'radius': my_attrs['R_p (Rearth)'] = items[4]
                            elif items[2] == 'mass': my_attrs['M_p (Mearth)'] = items[4]
                        elif items[1] == 'Equatorial': my_attrs['Wind speed kms'] = items[8]
                        elif items[1] == 'Rotation': my_attrs['Rotation period sec'] = items[4]
                        elif items[1] == 'Number': my_attrs['Number of rotations'] = items[6]
                        elif items[1] == 'Diffusion': my_attrs['vertical mixing'] = items[5]
                        elif items[1] == 'Photochemistry': my_attrs['photochemistry'] = items[5]
                else:
                    break
        return my_attrs   


    def read_2D_pac(self,
                    multiplexing = True,
                    xarray_format = True,
                    only_last_rot=False):
        
        """Function that reads in the output of the PAC2D simulation. This is the .dat file
        that is produced in the output directory of the PAC2D simulation. This function returns
        a xarray Dataset object, containing all the variables in the .dat file. The dimensions of
        the xarray Dataset are: (rotation, longitude, pressure)

        Args:
            only_last_rot (bool, optional): Boolean that allows you to select only the last rotation.
                                            Defaults to False.

        Raises:
            Exception: If there is more than one .dat file in the output directory, this function
                        will raise an exception. This is because there should only be one .dat file
                        in the output directory. If there are more, this is probably a mistake.
        """

        #Get the paths of all the _out_*.dat files in the PAC directory
        paths_dat_files = [self.outputdir + file for file in os.listdir(self.outputdir) if (file.endswith('.dat') and "_out_" in file )]
        if len(paths_dat_files) != 1:
            raise Exception(f"More than one .dat file found in {self.outputdir}. Double check this... For now, we continue with the first file.")
        f0 = paths_dat_files[0] #There should only be one output file in the pac directory...

        #Read in the pac data file
        df_pac = load_the_PAC2D_dat_file(f0)
        #Read in the attributes of the pac data file
        my_attrs = self._read_my_attrs(f0)

        # Change the name of the [] coordinates to _ because 
        # this is more easy for python to interpret and call in the xarray framework
        df_pac.rename(columns={
            'longitude[pi]':'longitude_pi',
            'height[km]':'height_km',
            'pressure[bar]':'pressure_bar',
            'temperature[K]':'temperature_K'
            },inplace=True
        )

        #Define the dimensions of our dataset
        nb_rotations = np.int(np.max(df_pac['longitude_pi']) // 2) #2pi means 1 rotations

        #Throw away the very last longitude because this is actually the start of rotation: nb_rotations+1
        #for example, 3 rotations of simulation would result in the very last longitude to be called rotation 4
        #Since the code converges to a steady state, this should not matter of you would include [0,2pi[ or ]0,2pi]
        #We choose [0, 2pi] because this is the most intuitive way of thinking about it.
        df_pac = df_pac[df_pac["longitude_pi"] < nb_rotations*2]

        #Add column that maps longitudes from 0 to 2(pi) with nb_vertical repititions of the longitude
        longs_0_2pi = sorted(df_pac['longitude_pi'][df_pac['longitude_pi'] < 2 ])

        #Check if we only keep have to keep the last rotation
        if only_last_rot:
            df_pac = df_pac[df_pac["longitude_pi"] >= (nb_rotations-1)*2]
            nb_rotations = 1

        #Repeat the list nb_rotations times to fill the column
        longs_0_2pi = longs_0_2pi * nb_rotations
        if 'longitude_npi' not in df_pac.columns:
            df_pac.insert(loc=0, column='longitude_npi', value=longs_0_2pi)

        #Add column that maps the rotations from 0 to nb_rotations with nb_rotations repititions of the rotations number
        rot_nbs = np.int0(df_pac['longitude_pi'] // 2)
        if 'rotation' not in df_pac.columns:
            df_pac.insert(loc=0, column='rotation', value=rot_nbs)

        ## APPLY MULTIPLEXING (I forgot what this does)
        if multiplexing:
            df_multi = df_pac.set_index(['rotation', 'longitude_npi', 'pressure_bar'])
        else:
            df_multi = df_pac

        if xarray_format:
            #Set the xarray dataset
            self.pac2d_xr = xr.Dataset.from_dataframe(df_multi)
            self.resolution_z = self.pac2d_xr.dims["pressure_bar"]
            self.resolution_l = self.pac2d_xr.dims["longitude_npi"]
        else:
            self.pac2d_xr = df_multi
            self.resolution_z = len(set(df_multi["pressure_bar"]))
            self.resolution_l = len(set(df_multi["longitude_npi"]))
            # self.resolution_z = len(df_multi.index.get_level_values("pressure_bar").unique())
            # self.resolution_l = len(df_multi.index.get_level_values("longitude_npi").unique())

        #Set the xarray dataset
        # self.pac2d_xr = xr.Dataset.from_dataframe(df_multi)

        #Now also read a bunch of metadata
        self.my_attrs = my_attrs
        self.R_pl = float(self.my_attrs["R_p (Rearth)"]) * cte.R_earth.cgs.value
        self.M_pl = float(self.my_attrs["M_p (Mearth)"]) * cte.M_earth.cgs.value
        self.R_star = float(self.my_attrs["R_star (Rsun)"]) * cte.R_sun.cgs.value
        self.gravity = self._compute_gravity()
        self.rotation_period = float(self.my_attrs["Rotation period sec"])

    def get_pac2d_xr(self):
        return self.pac2d_xr
    
    def get_rotation_period(self):
        return self.rotation_period

    def get_subset_xrdata(self, side = "substellar", width_pi = 0.1, rotation = "last"):
        """This function returns a subset of the xarray Dataset, containing only the data
        for a specific side of the planet. The side can be selected by the center of the side
        in longitude and the width of the side in pi. For example, the substellar point is
        located at 0.0 longitude and could have a width of 0.1 pi. This means that the substellar
        point is located at 0.0 longitude and the evening point is located at 0.5 longitude.

        Args:
            side (str, optional): String that describes the longitude in words . Defaults to "substellar".
            width_pi (float, optional): Width of the interval around the selected longitude in fractions of pi. Defaults to 0.1.
            rotation (str, optional): Does not have to be a string. Defaults to "last".
        """        

        #Select rotation number if the given rotation variable is a string
        if rotation == "last":
            rotation = max(self.pac2d_xr.rotation)

        #Select the rotation
        if rotation!=None:
            xr_2D = self.pac2d_xr.sel(rotation = rotation)
        else:
            xr_2D = self.pac2d_xr
            print("[-I] No rotation was selected, so all rotations are included in the subset")

       #Select side in longitude
        if (type(side) != str):
            c = side
        elif side == "substellar":
            c = 0
        elif side == "evening":
            c = 0.5
        elif side == "antistellar":
            c = 1
        elif side == "morning":
            c = 1.5
        elif side == "dayside":
            c = 0
            width_pi = 0.99
            print("[-I] Dayside is selected, so the width is set to 0.99")
        elif side == "nightside":
            c =1
            width_pi = 0.99
            print("[-I] Nightside is selected, so the width is set to 0.99")
        else:
            print(f"[-I] Side {side} not recognised. Returning full dataset.")
            return xr_2D

        return get_subset_around_long(xr_2D, center = c, width=width_pi)


    def get_list_of_molecules(self):
        """Function that returns a list of molecules that are present in the PAC2D simulation

        Returns:
            _type_: _description_
        """
        removables = ['longitude_pi', 'height_km', 'temperature_K']

        all_variables = list(self.pac2d_xr.keys())

        for removable in removables:
            try:
                all_variables.remove(removable)
            except:
                print(f"Unable to remove {removable} from list of variables")

        return all_variables


    def get_mean_molecular_weight(self, side ="substellar", width_pi=0.1, rotation='last', per_layer = False):
        """Function that computes the mean molecular weight of the atmosphere. 

        Args:
            side (str, optional): Specify which longitude you want to use. Defaults to "substellar".
            width_pi (float, optional): Width of the slice around the specified longitude. Defaults to 0.1.
            rot (str, optional): Rotation number to include in the slice. Defaults to 'last'.
            per_layer (bool, optional): Boolean that controls if the mean molecular weight is computed per layer or not. Defaults to False.
        """        

        #Number fractions: Take the mean over the longitude and altitude for the last rotation
        # nf_log = self.pac2d_xr.sel(rotation = max(self.PACdat.rotation))
    #     print(nf_log['H'].values)
        #Extract a list of molecules from the pac directory
        molecule_list = self.get_list_of_molecules()

        # nf = self.get_number_fractions(side=side, width_pi=width_pi, rot=rot)
        pac2d_xr_subset = self.get_subset_xrdata(side = side,
                                                  width_pi = width_pi,
                                                    rotation = rotation)


        total= 0.00
        for mol in molecule_list:

            #Get the weight of this molecule in amu
            mol_weight = get_molecular_weight_molecule(mol)

            #get the mean number fraction in the atmosphere
            nf_mol = np.power(10, pac2d_xr_subset[mol])

            #Choose the required output...
            if not per_layer:
                nf_mol = nf_mol.mean() #take mean over all dimensions per variable

            #transform to numpy floats
            nf_mol =  np.float64(nf_mol.values)

            total += mol_weight * nf_mol

        return total
    
    def get_lpt_file(self):
        return self.lpt_xr

    def load_lpt_file(self, keyword=''):
        """Function that loads the .lpt file that is used to init the PAC2D simulation.

        Args:
            keyword (str, optional): Keyword that is used to select the correct .lpt file. Defaults to ''.
        """

        path_file, _ = get_file_path_and_name(self.outputdir, ".lpt", keyword=keyword)

        with open(path_file, "r") as f:
            lines = f.readlines()
            pressures = np.fromstring(lines[4], dtype=np.float64, sep =' ')

        d =  np.loadtxt(path_file, skiprows=6)

        #longitudes
        longitudes_i = d[:,1]
        N = len(longitudes_i)
        longitudes_pi  = np.linspace(0, 2, N+1)[:-1]

        #temperature data
        data = d[:, 1:]

        #make a xarray object
        self.lpt_xr = xr.DataArray(data = data.T,
                             coords = {"pressure_bar": pressures, "longitude_npi": longitudes_pi},
                             dims = ["pressure_bar", "longitude_npi"],
                             name = 'temperature_K')

        self._compute_MWS()


    def _compute_MWS(self):
        """Computes the mean wind speed of the atmosphere
        """
        #Open lpt file
        lpt_file = [self.outputdir + pa for pa in os.listdir(self.outputdir) if pa.endswith('.lpt')][0]
        with open(lpt_file, "r") as f:
            mws = float(f.readline().split()[0])
            
        self.mws = mws

    def plot_LPT(self, degrees = False):
        """Plot the LPT file that shows the temperature as a function of pressure and longitude.


        Args:
            degrees (bool, optional): Boolean that controls the x-axis of the plot. Defaults to False.
        """        
        fig, ax = plt.subplots(1,1)
        lvls = 10

        xr_T  = self.pac2d_xr.sel(rotation = max(self.pac2d_xr.rotation))["temperature_K"]
        
        #Perform a roll to go from 0-2pi to -180*-180*
        if degrees==True:
            xr_T = roll_to_degrees(xr_T)
            
        #Determine the range of the colorscale
        vmax = np.max(xr_T.sel(longitude_npi = 0, pressure_bar = slice(1e-6,10)))
        
        #Plot the contourfill colors
        cnt = xr.plot.contourf(xr_T, ax=ax,
                            cmap=plt.cm.inferno,
                            x="longitude_npi",
                            y = "pressure_bar",
                            xincrease = True,
                            yincrease = False,
                            yscale = 'log',
                            levels = lvls*5,
                            robust = True,
                            vmax = vmax,
                            cbar_kwargs = {"label":r"Temperature [K]"}
                            )
        #Plot the lines on top of the contourfill
        im = xr.plot.contour(xr_T, ax=ax,
                            x="longitude_npi",
                            y = "pressure_bar",
                            xincrease = True,
                            yincrease = False,
                            yscale = 'log',
                            levels = lvls,
                            colors = 'gray',
                            robust = True,
                            vmax = vmax
                            )

        #General stuff to make the plot look nice
        ax.clabel(im, fmt='%2.0f', colors='white')
        ax.set_ylabel(r"Pressure [bar]")
        # ax.set_ylim([50, 5e-8])
        ax.set_title("")


        #Determine wether to plot in degrees or pi
        if degrees==True:
            ax.set_xlabel(r"Degrees [°]")
        else:
            ax.set_xlabel(r"Longitude [$\pi$]")
            

    def plot_chemistry_on_all_longitudes(self,
                                        list_molecules = ['NH3','HCN','CO2','H2O','CH4','CO','H2'],
                                        rotation = "last",
                                        ACE = None
                                        ):
        """Function that plots the chemistry of the atmosphere on all longitudes. The user can select
        which molecules to plot. The function also allows to plot the chemical equilibrium from the ACE

        Args:
            list_molecules (list, optional): List of molecules in PAC format to plot. Defaults to ['NH3','HCN','CO2','H2O','CH4','CO','H2'].
            rotation (str, optional): Rotation to use for the plot. Defaults to "last".
            ACE (_type_, optional): apu.ACE object holding the chemical equilibrium abundances. Defaults to None.
        """
        
        #Check if ACE is supplied and warn the user about the consequences
        if ACE != None:
            print("[-W] Plotting ACE chemical equilibrium. The user is responsible to know which longitudes the ACE file belongs to")

        #Select the rotation you want to plot
        if rotation =="last":
            rotation = np.max(self.pac2d_xr.rotation)
        xr_arr = self.pac2d_xr.sel(rotation = rotation)

        #Initialize the figure
        fig, ax = plt.subplots(1,1)
        cmap = plt.get_cmap("tab20")
        colors = list(cmap(np.linspace(0.01,1, len(list_molecules))))
        c_iter = iter(colors)
        
        #Plot the data
        cnt = 0 #counter to keep track of the longitude
        labelflag=True
        for long in xr_arr.longitude_npi:
            # if cnt%1==0:
            xr_long = xr_arr.sel(longitude_npi = long)
            for molecule in list_molecules:

                #Don't duplicate labels
                label = None
                if labelflag:
                    label = ppmol(molecule)

                #Plot the data
                c = next(c_iter)
                xr.plot.line(xr_long[molecule], ax = ax, y = 'pressure_bar', lw=0.5,
                            yscale = 'log', yincrease = False, label=label, c = c)
                
                #Plot the chemical equilibrium from ACE
                if ACE != None:
                    xr.plot.line(ACE.get_ace_xr()[molecule],
                                ax = ax,
                                y = 'pressure_bar',
                                yscale = 'log',
                                yincrease = False,
                                c = c,
                                linestyle = "--", # (0, (1,10)),
                                alpha = 0.6
                                )
                
            #set the labels to None from here onwards and reset the color iterator    
            labelflag=False
            c_iter = iter(colors)
            # cnt += 1        

        ax.grid(False)
        ax.set_xlabel("Molar fraction")
        ax.set_ylabel(r"Pressure [bar]")
        ax.set_xlim([-15, 1])
        legend = ax.legend(loc="lower left", labelcolor='linecolor',
                        handletextpad=0, handlelength=0, 
                        framealpha = 0.5, frameon=False)#bbox_to_anchor=(1.25, 1))
        for item in legend.legendHandles:
            item.set_visible(False)
    
        ax.set_title("")


    def plot_chemistry_on_4_longitudes(self,
            list_molecules = ['NH3', 'HCN', 'CO2', 'H2O', 'CH4', 'CO'],
            rotation = "last",
            width_pi = 0.01,
            ACE_dict = None
            ):
        """Plot the chemistry of the atmosphere on 4 longitudes. The user can select
        which molecules to plot. The function also allows to plot the chemical equilibrium from the ACE

        Args:
            list_molecules (list, optional): List of molecules in PAC format to plot. Defaults to ['NH3', 'HCN', 'CO2', 'H2O', 'CH4', 'CO'].
            rotation (str, optional): Rotation to use for the plot. Defaults to "last".
            width_pi (float, optional): Width of the column around each of the four longitudes. Defaults to 0.01.
            ACE_dict (_type_, optional): Dictionary containing the ACE objects for each of the four longitudes. Defaults to None.
        """


        sides = ['substellar', 'evening', 'antistellar', 'morning']
        
        #Check if ACE is supplied and warn the user about the consequences
        if ACE_dict != None:
            print("[-W] Plotting ACE chemical equilibrium. The user is responsible to know which longitudes the ACE file belongs to")

        #Select the rotation you want to plot
        # if rotation =="last":
        #     rotation = np.max(self.pac2d_xr.rotation)
        # xr_arr = self.pac2d_xr.sel(rotation = rotation)
        
        #Initialize the figure and colors
        fig, axs = plt.subplots(1,len(sides),
                                figsize = (3*len(sides), 4),
                                sharey=True, gridspec_kw = {'wspace':0.05}
                                )
        
        axs = axs.flatten()
        axs_iter = iter(axs)
        cmap = plt.get_cmap("tab20")
        colors = list(cmap(np.linspace(0.01,1, len(list_molecules))))       
        c_iter = iter(colors)
        
        #loop over the pac outputs        
        for j in range(len(sides)):
            side = sides[j]
            ax = next(axs_iter) 
            
            subset_xr = self.get_subset_xrdata(side = side,
                                                        width_pi = width_pi,
                                                        rotation = rotation)

            for molecule in list_molecules:
                c = next(c_iter)
                label = ppmol(molecule)

                # line = subset.sel(rotation = rot)[molecule]
                xr.plot.line(subset_xr[molecule], ax = ax,
                            y = 'pressure_bar',
                            yscale = 'log',
                            yincrease = False,
                            c = c,
                            label = label,
                            lw=2)
                
                #Plot the chemical equilibrium from ACE
                if side in ACE_dict:
                    if ACE_dict[side] != None:
                        ACE = ACE_dict[side]
                        xr.plot.line(ACE.get_ace_xr()[molecule], ax = ax,
                                    y = 'pressure_bar',
                                    yscale = 'log', yincrease = False,
                                    c = c, linestyle = "--", alpha = 1, lw=2)

            c_iter = iter(colors)
            ax.set_xlabel("Log molar fraction")
            if j==0:
                ax.set_ylabel(r"Pressure [bar]")
            else:
                ax.set_ylabel("")                
            ax.set_title("")
            ax.text(0.05, 0.94,
                f"{ppside(side)}",
                horizontalalignment="left",
                verticalalignment="bottom",
                clip_on=True,
                transform = ax.transAxes)

            # ax.set_ylim([50, 5e-8])
            ax.set_xlim([-15, 1])
            # ax.set_xticks([-14, -10, -6, -2])
            
        legend = axs[0].legend(loc="lower left", labelcolor='linecolor',
                        handletextpad=0, handlelength=0, 
                        # framealpha = 0.5,
                          frameon=False,
                                #    bbox_to_anchor=(1.4, 0.25)
                            )
        for item in legend.legendHandles:
            item.set_visible(False)


#            ____  ______
#     ____  / __ \/_  __/
#    / __ \/ /_/ / / /   
#   / /_/ / _, _/ / /    
#  / .___/_/ |_| /_/     
# /_/                    for petitRADTRANS usage and mass fractions


    def set_mass_fractions_and_PT(self, 
                           PRT_ATM,
                           additional_massFractions=None,
                           side = "evening",
                           width_pi=0.01, 
                           rotation = "last"):
        """function that computes the mass fractions of the atmosphere. This function
        also computes the mean molecular weight of the atmosphere and stores the temperature
        and pressure in a variable.

        Args:
            PRT_ATM (apu.PRT_ATM object): apu.PRT_ATM object that contains the petitRADTRANS atmosphere object.
            additional_massFractions (_type_, optional): Additional mass fractions that are not included in the PAC2D simulation. Defaults to None.
            side (str, optional): Longitude of the simulation. Defaults to "evening".
            width_pi (float, optional): Width of the slice around the selected longitude in fractions of pi. Defaults to 0.01.
            rotation (str, optional): Rotation number to include in the slice. Defaults to "last".

        Raises:
            Exception: If the mass fractions of the additional molecules are not included in the PAC2D simulation, this function will raise an exception.
        """
        
        #Mean molecular weight of the atmosphere
        mmw = self.get_mean_molecular_weight(side =side,
                                            width_pi=width_pi,
                                            rotation=rotation,
                                            per_layer = True)
        
        #Get subset of the xarray data at the correct longitude
        pac2d_xr_subset = self.get_subset_xrdata(side = side,
                                        width_pi = width_pi,
                                        rotation = rotation)
        
        mass_fractions = {}
        mass_fractions_missing = []

        #Species for which mass fractions should be computed
        mf_molecules = PRT_ATM.give_mf_molecules()

        for m in mf_molecules:

            try:
                
                #Apply formula of pRT to compute mass fractions from number fractions
                mu_i = get_molecular_weight_molecule(m)

                #Read the number fraction in log-format
                nf = pac2d_xr_subset[m].to_numpy()
                
                #Spit out the correct pRT name for the molecule e.g. SO2 (pac) becomes SO2_ExoAmes (pRT)
                pRT_key = translate_PAC_to_pRT_molecule(m)
                
                #Add to mass fraction dictionary
                mass_fractions[pRT_key] = (mu_i)/mmw * np.float64(np.power(10, nf))

            except:
                mass_fractions_missing.append(m)

        if additional_massFractions != None:
            for m in list(additional_massFractions.keys()):            
                mass_fractions[m] = additional_massFractions[m]
                mass_fractions_missing.remove(m)

        if len(mass_fractions_missing)!=0:
            raise Exception("The following molecules are not included in the mass fractions:", mass_fractions_missing)

        self.mass_fractions = mass_fractions
        self.PRT_ATM = PRT_ATM
        self.MMW = mmw
        self.temperatures = pac2d_xr_subset["temperature_K"].to_numpy()
        self.pressures = pac2d_xr_subset["pressure_bar"].to_numpy()



    def plot_convergence(self,
                        side = "substellar",
                        pressure = 1e-2,
                        molecule = 'NH3'):
        """Function that plots the convergence of the mixing ratio of a certain molecule as a function of the rotation number.
        This should be a straight line. If it is not, the simulation has not converged yet.

        Args:
            side (str, optional): Longitude of the simulation. Defaults to "substellar".
            pressure (_type_, optional): Pressure value at which to evaluate the convergence. Defaults to 1e-2.
            molecule (str, optional): Target molecule abundance. Defaults to 'NH3'.
        """        
        
        #compute the correct slicen of pressure and longitude
        xr_arr = self.get_subset_xrdata(side = side,
                                    width_pi = 1e-1,
                                    rotation = None)
        xr_arr = xr_arr.sel(pressure_bar=pressure, method='nearest')

        #Initialize the figure
        fig, ax = plt.subplots(1,1,figsize=(8,4),
                                sharey=False, sharex=True, 
                                gridspec_kw = {'hspace':0.02}
                                )

        #Remove the first datapoint because this is the initial guess
        x = xr_arr.rotation.values[1:]
        y = 10**(xr_arr[molecule].values[1:])
        ax.plot(x, y, marker="*")
        
        #Set the labels
        ax.set_xlabel("Rotation number")
        ax.set_ylabel("Log molar fraction")
        ax.ticklabel_format(style='sci', axis='y',
                            useMathText=True, scilimits=[2,3], useOffset=False)
        # ax.legend(loc="best")
        ax.grid(True, alpha=0.2)            


#------------------------------------------------------------------------------------------------------
# other functions
#------------------------------------------------------------------------------------------------------

def load_the_PAC2D_dat_file(file_path):
    """Function that loads the PAC2D .dat file and returns a pandas dataframe

    Args:
        file_path (_type_): Path to the .dat file

    Returns:
        _type_: Dataframe containing the data from the .dat file
    """
    f = file_path

    #Get the top and end indices for the rows that start with a '!'
    top, end = give_skiprow_top_and_end(f)

    #First extract a list of row that start with "!", then select the last row [-1] and cut of the "! " by [2:]
    headers_string = [line for i, line in enumerate(open(f)) if line.startswith('!')][-1][2:]

    # Split the different variable names in headers_string
    headers_list = headers_string.split()

    #Keep none if there are no final rows that start with '!'
    max_row = None
    if len(end) != 0:
        max_row = end[0] - len(top)

    #Extract the data without column names
    data = np.loadtxt(f, skiprows= len(top), max_rows=max_row)

    # #Include or exclude the altitude, pressure, temperature headers (apt) accoording to prefference
    # if not include_apt_headers:
    #     headers_list.remove('longitude[pi]')
    #     headers_list.remove('height[km]')
    #     headers_list.remove('pressure[bar]')
    #     headers_list.remove('temperature[K]')
    #     data = data[:,3:]

    #Make a pandas dataframe
    df = pd.DataFrame(data = data,
                     columns= headers_list,
                     dtype= np.float64)

    return df