import pandas as pd
import re
import numpy as np
from matplotlib import pyplot as plt 
from scipy.optimize import curve_fit
from matplotlib.ticker import ScalarFormatter
from datetime import datetime
import cmath

from datetime import datetime
import numpy as np
import cmath

class ModelManager:
    def __init__(self, model):
        """
        Initialize the ModelManager class with a COMSOL model instance.
        
        Parameters:
        model (object): The COMSOL model object that the class will manage.
        """
        self.model = model

    def run(self, study):

        ida = self.model.parameter('ida')
        coren = self.model.parameter('nCore', evaluate=1)
        cladn = self.model.parameter('nClad', evaluate=1)
        pml = self.model.parameter('pml_extend')
        air_offset = self.model.parameter('air_offset')
        idaPLM = self.model.parameter('idaPLM', evaluate=1)
        rmodel = self.model.parameter('radius')

        print('Model return:\n')
        print(f'ida: {ida}')
        print(f'nClad: {cladn}')
        print(f'nCore: {coren}')
        print(f'pml_um: {pml}')
        print(f'air_offset: {air_offset}')
        print(f'idaPLM({rmodel}, {ida}): {idaPLM}')

        self.model.solve(study) 
        


    def set_wavelength(self, lambda_nm, cladding, core, pml_um=-1, air_offset_um=-1, note=None):
        """
        Set the wavelength and material properties for the COMSOL model.
        
        Parameters:
        lambda_nm (float): The wavelength in nanometers.
        cladding (str): Path to the file containing the cladding material data.
        core (str): Path to the file containing the core material data.
        pml_um (float): The PML extension in micrometers (default: -1).
        air_offset_um (float): The air offset in micrometers (default: -1).
        note (str): Optional note to include in the description.
        """
        # Set the time stamp and description
        if note is not None:
            self.model.description('readme', f'{note} @{lambda_nm}[nm]')
        else:
            note = datetime.now().strftime("%m/%d/%Y %I:%M%p")
            try:
                self.model.description('readme', f'{note} @{lambda_nm}[nm]')
            except ValueError:
                check = self.model.description('readme')
                print(f'readme has been set to: {check}\n')

        # Prepare the values
        nClad = self.n(lambda_nm, cladding)[1]
        nCore = self.n(lambda_nm, core)[1]

        # Load parameters into the COMSOL model
        self.model.parameter('nCore', f'{nCore}')
        self.model.parameter('nClad', f'{nClad}')
        self.model.parameter('ida', f'{lambda_nm}[nm]')

        if pml_um > 0:
            self.model.parameter('pml_extend', f'{pml_um}[um]')
        if air_offset_um > 0:
            self.model.parameter('air_offset', f'{air_offset_um}[um]')

        # Print the comparison between the set and retrieved parameters
        ida = self.model.parameter('ida')
        coren = self.model.parameter('nCore', evaluate=1)
        cladn = self.model.parameter('nClad', evaluate=1)
        pml = self.model.parameter('pml_extend')
        air_offset = self.model.parameter('air_offset')
        idaPLM = self.model.parameter('idaPLM', evaluate=1)
        rmodel = self.model.parameter('radius')

        print('You want:\n')
        print(f'ida: {lambda_nm}')
        print(f'nClad: {nClad}')
        print(f'nCore: {nCore}')

        print('Model return:\n')
        print(f'ida: {ida}')
        print(f'nClad: {cladn}')
        print(f'nCore: {coren}')
        print(f'pml_um: {pml}')
        print(f'air_offset: {air_offset}')
        print(f'idaPLM({rmodel}, {ida}): {idaPLM}')

    def export(self, table_name, export_path):
        """
        Export a table with metadata and headers from the COMSOL model.

        Parameters:
        table_name (str): The name of the table to export.
        export_path (str): The file path where the exported file will be saved.
        """
        # Step 1: Access the table node dynamically by name
        result = self.model/'tables'/table_name
        
        # Step 2: Save the table to a file
        result.java.saveFile(export_path)
        
        # Step 3: Retrieve the headers and convert them to Python strings
        headers = [str(header) for header in result.java.getColumnHeaders()]
        
        # Step 4: Dynamically retrieve the model name and table name
        model_name = str(self.model.name())  # Fetching model name dynamically
        table_display_name = str(result.name())  # Fetching table name dynamically

        # Step 5: Open the saved file, add metadata, headers with '%' and save again
        with open(export_path, 'r') as file:
            content = file.readlines()

        # Step 6: Prepare the metadata and headers to be added
        metadata = (
            f"% Model:              {model_name}.mph\n"
            f"% Date:               {datetime.now().strftime('%b %d %Y, %I:%M%p')}\n"
            f"% Table:              {table_display_name}\n"
        )

        # Add '%' to the start of the headers line
        headers_line = '%' + '\t'.join(headers) + '\n'

        # Step 7: Combine metadata, headers, and file content
        content.insert(0, metadata)  # Insert metadata at the top
        content.insert(1, headers_line)  # Insert headers after metadata

        # Step 8: Write the modified content back to the file
        with open(export_path, 'w') as file:
            file.writelines(content)

        print(f"Table saved with metadata and headers at: {export_path}")

    def n(self, lambda_nm, material_file):
        """
        Load refractive index data from a material file and return the closest index to lambda_nm.

        Parameters:
        lambda_nm (float): The wavelength in nanometers.
        material_file (str): Path to the file containing the material data.

        Returns:
        Tuple: Closest wavelength and corresponding refractive index
        """
        data = np.loadtxt(material_file)
        # If the wavelength data is already in nanometers, don't convert it
        if data[:, 0][0] > 100:
            wavelength_nm = data[:, 0]
        else:
            wavelength_nm = data[:, 0] * 1e9  # Convert to nanometers if in meters

        refractive_index = data[:, 1]
        closest_idx = np.argmin(np.abs(wavelength_nm - lambda_nm))
        return wavelength_nm[closest_idx], refractive_index[closest_idx]

    def beta_r(self, radius_um, offset_um, lambda_nm, core, cladding, width_core_um=0):
        """
        Compute beta_r for the given parameters.
        
        Parameters:
        radius_um (float): Radius in micrometers.
        offset_um (float): Offset in micrometers.
        lambda_nm (float): Wavelength in nanometers.
        core (str): Path to the file containing the core material data.
        cladding (str): Path to the file containing the cladding material data.
        width_core_um (float): Width of the core in micrometers (default: 0).
        
        Prints the result.
        """
        ida = lambda_nm * 1e-9  # Convert to meters
        nClad = self.n(lambda_nm, cladding)[1]  # Cladding refractive index
        nCore = self.n(lambda_nm, core)[1]  # Core refractive index
        dr = (offset_um + width_core_um / 3) * 1e-6  # Convert to meters
        radius = radius_um * 1e-6  # Convert to meters

        # Calculate the expression
        inside_sqrt = nClad**2 - (nCore * radius / (radius + dr))**2
        result = ida / cmath.sqrt(inside_sqrt)  # Use cmath.sqrt to handle complex numbers
        result_um = result * 1e6  # Convert back to micrometers

        minairoff = (nCore / nClad - 1) * radius_um - width_core_um / 3
        print(f"Ida_Radial({lambda_nm}) = {result_um} [um]")
        print(f"minimum air offset is {minairoff:.4f} [um]")

    """
    Export a table with metadata, headers, and properly formatted output from COMSOL using Python mph.

    Parameters:
    model (object): The COMSOL model object.
    table_name (str): The name of the table to be exported.
    export_path (str): The path where the exported file will be saved.
    """
    # Step 1: Access the table node dynamically by name
    result = model/'tables'/table_name
    
    # Step 2: Save the table to a file
    result.java.saveFile(export_path)
    
    # Step 3: Retrieve the headers and convert them to Python strings
    headers = [str(header) for header in result.java.getColumnHeaders()]
    
    # Step 4: Dynamically retrieve the model name and table name
    model_name = str(model.name())  # Fetching model name dynamically
    table_display_name = str(result.name())  # Fetching table name dynamically

    # Step 5: Open the saved file, add metadata, headers with '%' and save again
    with open(export_path, 'r') as file:
        content = file.readlines()

    # Step 6: Prepare the metadata and headers to be added
    metadata = (
        f"% Model:              {model_name}.mph\n"  # Dynamically retrieved model name
        f"% Date:               {datetime.now().strftime('%b %d %Y, %I:%M%p')}\n"
        f"% Table:              {table_display_name}\n"  # Dynamically retrieved table name
    )

    # Add '%' to the start of the headers line
    headers_line = '%' + '\t'.join(headers) + '\n'

    # Step 7: Combine metadata, headers, and file content
    content.insert(0, metadata)  # Insert metadata at the top
    content.insert(1, headers_line)  # Insert headers after metadata

    # Step 8: Write the modified content back to the file
    with open(export_path, 'w') as file:
        file.writelines(content)

    print(f"Table saved with metadata and headers at: {export_path}")