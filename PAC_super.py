from .utils import *
from .petitRADTRANS_init import *
import astropy.constants as cte
import matplotlib.pyplot as plt

class PAC_super:

    def __init__(self, outputdir, path_figsave=None):
        """Class that is the super class of the PAC1D and PAC2D classes. 
        It contains the basic functions that are used in both classes.

        Args:
            outputdir (str): Path to the directory that contains the input and output files of the PAC simulation
            path_figsave (str, optional): Path to the directory where the figures need to be saved. Defaults to None.
        """

        #Set the plotting directory
        self.path_figsave = path_figsave        
        if path_figsave != None:
            if not os.path.isdir(path_figsave):
                self._create_folder()

        #Variables for the output of the pac2D simulation
        self.outputdir = outputdir

        #Data of the pac directory
        self.my_attrs = None
        self.R_pl = None              #Radius in R_earth
        self.M_pl = None              #Mass in M_earth
        self.R_star = None            #Radius of the star in R_sun
        self.gravity = None           #Gravity in cgs

        #Parameters for the transmission spectrum with petitRADTRANS
        self.mass_fractions = None
        self.temperatures = None      #Temperature in K
        self.pressures = None         #Pressure in bar   
        self.MMW = None               #Mean molecular weight in amu
        self.PRT_ATM = None
        self.x = None                 #Wavelength in um
        self.y = None                 #Transit depth in %

        # Constructor code here
        self.star_df = None
        self.eddy_df = None

    def _create_folder(self):
        """Function that creates a folder to save the figures in and checks if the folder already exists

        Raises:
            ValueError: Error that is raised when the non-existing folder cannot be created
        """
        try:
            os.makedirs(self.path_figsave)
            print(f"Folder '{self.path_figsave}' created successfully.")
        except OSError as e:
            raise ValueError(f"Unable to create folder {self.path_figsave}. Error: {e}")

    def get_outputdir(self):
        return self.outputdir
    
    def get_my_attrs(self):
        return self.my_attrs
    
    def get_R_pl(self):
        return self.R_pl

    def get_m_PL(self):
        return self.M_pl
    
    def get_R_star(self):
        return self.R_star

    def _compute_gravity(self):
        """Function that computes the gravity (cgs) of the planet for which the simulation is performed

        Returns:
            float: Gravity in cgs-units of the hypothetical planet, for which the ace/pac simulation is perfomed.
        """
        R_pl = self.R_pl #* cte.R_earth.cgs.value
        M_pl = self.M_pl #* cte.M_earth.cgs.value
        G = cte.G.cgs.value

        gravity = (G * M_pl)/(R_pl**2)

        return gravity
    
    def get_gravity(self):
        return self.gravity
    

#            ____  ______
#     ____  / __ \/_  __/
#    / __ \/ /_/ / / /   
#   / /_/ / _, _/ / /    
#  / .___/_/ |_| /_/     
# /_/                    petitRADTRANS usage and mass fractions

    def set_mass_fraction_to_zero(self, mols=[]):
        """Set the mass fractions of certain molecules to zero

        Args:
            mols (list, optional): List of molecules which need to be set to zero. Defaults to [].
        """        
        for m in mols:

            m_pRT = translate_PAC_to_pRT_molecule(m)
            self.mass_fractions[m_pRT] = np.zeros_like(self.mass_fractions[m_pRT])

    def keep_mass_fraction_only_of(self, mols=[]):
        """Keep the mass fractions of certain molecules, set the others to zero

        Args:
            mols (list, optional): List of molecules for which we need to keep the mass fractions. Defaults to [].
        """        
        new_mass_fracs = {}
        mf_keys = list(self.mass_fractions.keys())
        for m in mf_keys:
            if m in mols:
                new_mass_fracs[m] = self.mass_fractions[m]
            else:
                new_mass_fracs[m] = np.zeros_like(self.mass_fractions[m])


    def compute_transmission_spectrum(self, 
                                      ref_R = 1,
                                      radius={}, sigma_lnorm=1.05,
                                      Pcloud=None, haze_factor=None,
                                      contribution=False):
        """Function that computes the transmission spectrum of the planet for which the PAC simulation is performed

        Args:
            ref_R (int, optional): Radius at which the amospheric reference pressure 0.01 bar is reached . Defaults to 1.
            radius (dict, optional): Dictionary containing the radius distribution of the cloud particles in function of the pressure. Defaults to {}.
            sigma_lnorm (float, optional): Distribution of cloud particle size. Defaults to 1.05.
            Pcloud (float, optional): Base pressure of cloud deck. Defaults to None.
            haze_factor (float, optional): Haze factor. Defaults to None.
            contribution (bool, optional): Flag that allows the user to compute the contribution function the transmission spectrum. Defaults to False.
        """        


        self.PRT_ATM.set_up(self.pressures)
    
        self.PRT_ATM.atmosphere.calc_transm(

            #Temperature structure
            temp = self.temperatures,

            #MAss fractions dictionary
            abunds = self.mass_fractions,

            #Surface gravity
            gravity = self.gravity,

            #mean molecular weight
            mmw = self.MMW,

            #Reference radius in cm
            R_pl = ref_R * self.R_pl,

            #Reference pressure where R(P=P0) = R_pl
            P0_bar=0.01,
            
            #clouds particle distribution
            radius = radius,
            sigma_lnorm = sigma_lnorm,
            
            #Simple cloud stuff
            haze_factor = haze_factor, #10,
            
            Pcloud = Pcloud, #0.01
            
            contribution=contribution
            
        )

        c = cte.c.cgs.value #cm/s
        x = c/self.PRT_ATM.atmosphere.freq/1e-4                               #Wavelength in um
        y = (self.PRT_ATM.atmosphere.transm_rad/self.R_star)**2 *100          #Transit depth in %

        print("[-I] Succesfully computed the transmission spectrum")

        self.x = x
        self.y = y

    def get_transmission_spectrum(self):
        return self.x, self.y

    def plot_transmission_spectrum(self, ax=None, log=False):
        """Function that plots the transmission spectrum

        Args:
            ax (matplotlib.ax object, optional): axis to plot the transmission spectrum on. Defaults to None.
            log (bool, optional): Boolan to choose linear/log scale of x-axis. Defaults to False.
        """

        if ax==None:
            fig, ax = plt.subplots()
            
        ax.plot(self.x, self.y, label="PRT")

        if log:
            ax.set_xscale("log")
        ax.set_xlabel(r"Wavelength ($\mu$m)")
        ax.set_ylabel(r'Transit depth ($\%$)')


    #Stellar spectrum
    def load_star_file(self, keyword=''):
        """Function that loads in the stellar spectrum

        Args:
            keyword (str, optional): Keyword that should be in the name of the .star file. Defaults to ''.
        """
        self.star_df = load_file_to_dataframe(self.outputdir , ext = ".star", keyword=keyword)

    def get_star_file(self):
        return self.star_df
    
    def plot_star_spectrum(self, ax=None, log=True):
        """Function that plots the stellar spectrum

        Args:
            ax (matplotlib.ax object, optional): Axis to plot the stellar spectrum on. Defaults to None.
            log (bool, optional): _description_. Defaults to True.
        """

        headers_list = list(self.star_df.columns)
        if ax==None:
            fig, ax = plt.subplots()
            
        ax.semilogy(self.star_df[headers_list[0]],
                self.star_df[headers_list[1]],
                  label="Star")

        if log:
            ax.set_xscale("log")

        ax.set_xlabel(r"Wavelength (nm)")
        ax.set_ylabel(r'Flux (erg s-1 cm-2 Hz-1 sr-1)')
    
    #Eddy diffusion coefficient
    def load_eddy_file(self, keyword=''):
        """Function that loads in the eddy diffusion coefficient profile

        Args:
            keyword (str, optional): Keyword that should be in the name of the .eddy file. Defaults to ''.
        """
        self.eddy_df = load_file_to_dataframe(self.outputdir , ext = ".eddy", keyword=keyword)
    
    def get_eddy_file(self):
        return self.eddy_df
    
    def plot_eddy_diffusion(self, ax=None, log=True):
        """Function that plots the eddy diffusion coefficient profile

        Args:
            ax (matplotlib.ax object, optional): Axis to plot the eddy diffusion coefficient profile on. Defaults to None.
            log (bool, optional): Boolean to choose linear/log scale of x-axis. Defaults to True.
        """

        headers_list = list(self.eddy_df.columns)

        if ax==None:
            fig, ax = plt.subplots()
            
        ax.semilogy(self.eddy_df[headers_list[1]],
                self.eddy_df[headers_list[0]],
                  label="Eddy diffusion")

        if log:
            ax.set_xscale("log")

        ax.set_ylabel(r"Pressure (bar)")
        ax.set_xlabel(r'Eddy diffusion coefficient (cm2 s-1)')
        ax.set_ylim(*ax.get_ylim()[::-1])
        ax.set_xlim([1e5, 1e13])