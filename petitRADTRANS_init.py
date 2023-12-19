#In the folder where you keep your petitRADTRANS opacities, go to 
# ./input_data/opacities/lines/corr_k/
#There you find the opacity files for the line species, which should connect the PAC names.

#Example: I downloaded an opacity file for SO2 from ExoMol, which is stored in the folder SO2_ExoAmes.
#         Therefore, I need to link the PAC name "SO2" to the pRT name "SO2_ExoAmes".

#Change this to match the correct files in your corr_k folder
### "PAC name" : "pRT name"

PAC_to_pRT = {
    ("SO", "SO"),
    ("C2H2", "C2H2"), 
    ("CO2", "CO2"), 
    ("OH", "OH"), 
    ("HCN", "HCN"), 
    ("NH3", "NH3"), 
    ("CH4", "CH4"), 
    ("CO", "CO_all_iso_HITEMP"), 
    ("H2O", "H2O_Exomol"), 
    ("H2S", "H2S"), 
    ("SO2", "SO2_ExoAmes")
}


def translate_PAC_to_pRT_molecule(PAC_name):
    """Function that translates the PAC name to the pRT name

    Args:
        PAC_name (str): Name of the molecule in the PAC format

    Returns:
        str: Name of the molecule in the pRT format
    """
    for PAC, pRT in PAC_to_pRT:
        if PAC == PAC_name:
            return pRT
    
    return PAC_name


def translate_pRT_to_PAC_molecule(pRT_name):
    """Function that translates the pRT name to the PAC name

    Args:
        pRT_name (str): Name of the molecule in the pRT format

    Returns:
        str: Name of the molecule in the PAC format
    """
    for PAC, pRT in PAC_to_pRT:
        if pRT == pRT_name:
            return PAC
    
    return pRT_name