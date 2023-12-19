from .utils import *
from .petitRADTRANS_init import *

import numpy as np
import scipy as sp
from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc
import matplotlib.pyplot as plt
import sys, os

import xarray as xr
import pandas as pd
import re

class PRT_ATM:
    def __init__(self, line_species = ['SO', 'C2H2', 'CO2', 'OH', 'HCN', 'NH3', 'CH4', 'CO', 'H2O', 'H2S', 'SO2'],
                       rayleigh_species = ["He", "H2"],
                       cloud_species = ['Mg2SiO4(c)_cd'],
                       continuum_opacities = ['H2-H2', 'H2-He'],
                       wlen_bords_micron = [0.9, 20]
                       ):
        """Class that makes a petitRADTRANS atmosphere object

        Args:
            line_species (list, optional): Lines species. Defaults to ['SO', 'C2H2', 'CO2', 'OH', 'HCN', 'NH3', 'CH4', 'CO', 'H2O', 'H2S', 'SO2'].
            rayleigh_species (list, optional): Rayleigh species. Defaults to ["He", "H2"].
            cloud_species (list, optional): Cloud species. Defaults to ['Mg2SiO4(c)_cd'].
            continuum_opacities (list, optional): Continuum opacities. Defaults to ['H2-H2', 'H2-He'].
            wlen_bords_micron (list, optional): Wavelength range in micron. Defaults to [0.9, 20].
        """
                       
        self.line_species = line_species
        self.rayleigh_species = rayleigh_species
        self.cloud_species = cloud_species
        self.continuum_opacities = continuum_opacities
        self.wlen_bords_micron = wlen_bords_micron
        self.atmosphere = None

    def get_line_species(self):
        return self.line_species
    
    def get_rayleigh_species(self):
        return self.rayleigh_species

    def get_cloud_species(self):
        return self.cloud_species

    def get_continuum_opacities(self):
        return self.continuum_opacities

    def get_wlen_bords_micron(self):
        return self.wlen_bords_micron

    def get_atmosphere(self):
        return self.atmosphere
    
    def load_atmosphere_object(self):
        """Function that loads in the atmosphere object
        """

        self.atmosphere = Radtrans(

        #line species: molecular or atomic line absorbers
        line_species = self.pac_to_prt(self.line_species),

        #rayleigh species: load the Rayleigh scattering cross sections
        rayleigh_species = self.rayleigh_species,

        #cloud_species: which cloud species to include`
        cloud_species = self.cloud_species,

        #continuum: load the collision induced absorption cross sections (CIA)
        continuum_opacities = self.continuum_opacities,

        #wavelength range in micron
        wlen_bords_micron = self.wlen_bords_micron
        
        )

    def pac_to_prt(self, mols):
        """Function that translates the PAC names to the pRT names

        Args:
            mols (list): List containing molecules in the PAC format

        Returns:
            list: List containing molecules in the pRT format
        """
        prt_format = []
        for mol in mols:
            prt_format.append(translate_PAC_to_pRT_molecule(mol))

        return prt_format

    def prt_to_pac(self, mols):
        """Function that translates the pRT names to the PAC names

        Args:
            mols (list): List containing molecules in the pRT format

        Returns:
            list: List containing molecules in the PAC format
        """
        pac_format = []
        for mol in mols:
            pac_format.append(translate_pRT_to_PAC_molecule(mol))

        return pac_format


    def set_up(self, pressures_np):
        """Function that sets up the atmosphere object with the correct pressures

        Args:
            pressures_np (numpy.array): Array of pressure values in bar in ascending order
        """
        try:
            self.atmosphere.setup_opa_structure(pressures_np)
        except AttributeError:
            print("[ERROR] You must first init a atmosphere object through the function `load_atmosphere_object` before you can run set_up()")

    def give_mf_molecules(self):
        """Function that returns a list of all molecules for which mass fractions should be computed

        Returns:
            list: All molecules for which mass fractions should be computed
        """

        mf_molecules = []

        mf_molecules.extend(self.line_species)     
        mf_molecules.extend(self.rayleigh_species)
        mf_molecules.extend(self.cloud_species)

        # print("The molecules for which mass fractions should be computed are:", mf_molecules)
        return mf_molecules