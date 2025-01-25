import pandas as pd
import re
import numpy as np
from matplotlib import pyplot as plt 
from scipy.optimize import curve_fit
from matplotlib.ticker import ScalarFormatter
from datetime import datetime
import cmath
from scipy.interpolate import griddata
import matplotlib.patches as patches
import plotly.graph_objects as go
import plotly.offline as pyo
import webbrowser
import os
def set_plot_style():
    """Configure global plot settings for a consistent and professional appearance."""
    plt.rcParams['figure.figsize'] = (8, 6)              # Set global figure size
    plt.rcParams['figure.dpi'] = 120                     # Set figure DPI for high resolution
    plt.rcParams['lines.linewidth'] = 2                  # Line width
    plt.rcParams['lines.markersize'] = 3                 # Marker size
    plt.rcParams['font.size'] = 12                       # Global font size
    plt.rcParams['legend.fontsize'] = 12                 # Legend font size
    plt.rcParams['legend.title_fontsize'] = 12           # Legend title font size
    plt.rcParams['axes.titlesize'] = 16                  # Title font size
    plt.rcParams['axes.labelsize'] = 14                  # Label font size
    plt.rcParams['grid.linestyle'] = '--'                # Gridline style
    plt.rcParams['grid.linewidth'] = 0.5                 # Gridline width
    plt.rcParams['axes.grid'] = True                     # Enable gridlines globally
    plt.rcParams['legend.loc'] = 'best'                  # Set best position for legend
def get_dog(positionr, position_neff, filepath, filter_ratio=False, ratio_col=None, threshold=0.9, extra_info_col=None, rank=0,ex_col=None,ez_col=None,vector_threshold=10,label_precision=1):
    radius_dict = {}
    extra_info_dict = {}  # Dictionary to store extra info associated with each radius

    # Read the file
    with open(filepath, 'r') as file:
        for line in file:
            # Skip lines that start with '%' or contain headers
            if line.startswith('%') or 'radius' in line or 'lambda' in line:
                continue

            # Split the line by whitespace
            data = line.split()

            # Ensure we have enough columns in the line to avoid index errors
            if len(data) > max(positionr, position_neff, ratio_col or 0, extra_info_col or 0):
                try:
                    # Check and apply the ratio filter if enabled
                    if filter_ratio and ratio_col is not None:
                        # Parse the ratio of core_energy / total_energy
                        ratio_value = float(data[ratio_col])
                        # Skip this row if the ratio is below the threshold
                        if ratio_value < threshold:
                            continue
                    if ex_col is not None:
                        ex_val = float(data[ex_col])
                        if ex_val<vector_threshold:
                            continue
                    if ez_col is not None:
                        ez_val = float(data[ez_col])
                        if ez_val<vector_threshold:
                            continue
                    # Extract the radius and real part of neff values
                    radius_value = round(float(data[positionr]), label_precision)  # Round radius to handle precision issues
                    neff_value = float(data[position_neff])  # Use only the real part of neff

                    # Optionally, extract the extra info if column is specified
                    extra_info = float(data[extra_info_col]) if extra_info_col is not None else None

                    # Check if the radius already exists in the dictionary
                    if radius_value in radius_dict:
                        # Append the neff value to the list for this radius
                        radius_dict[radius_value].append(neff_value)
                        if extra_info is not None:
                            extra_info_dict[radius_value].append(extra_info)
                    else:
                        # Initialize the entry with current neff in a list
                        radius_dict[radius_value] = [neff_value]
                        if extra_info is not None:
                            extra_info_dict[radius_value] = [extra_info]
                except ValueError:
                    # Skip lines where conversion fails
                    continue

    # Extract the desired ranked neff for each radius
    radii = []
    selected_neff = []
    selected_extra_info = []

    for radius, neff_list in radius_dict.items():
        # Sort neff values for each radius in descending order
        sorted_neff_list = sorted(neff_list, reverse=True)

        # Check if the requested rank exists
        if rank < len(sorted_neff_list):
            radii.append(radius)
            selected_neff.append(sorted_neff_list[rank])
            if extra_info_col is not None:
                selected_extra_info.append(extra_info_dict[radius][neff_list.index(sorted_neff_list[rank])])
        else:
            # Skip this radius if there are not enough neff values for the requested rank
            continue

    # Convert results to numpy arrays
    radii = np.array(radii)
    selected_neff = np.array(selected_neff)
    selected_extra_info = np.array(selected_extra_info) if extra_info_col is not None else None

    # Return radii, selected_neff, and optional extra_info
    return (radii, selected_neff, selected_extra_info) if extra_info_col is not None else (radii, selected_neff)
def write_array_with_name_to_txt(array, output_file_path, user_defined_name):
    """
    Writes a 1D NumPy array to a text file in the format: user_defined_name {value1 value2 ...}.
    
    Parameters:
    - array: 1D NumPy array to be written to the file.
    - output_file_path: String specifying the path to save the output file.
    - user_defined_name: String specifying the name to be used in the output file.
    """
    with open(output_file_path, 'w') as f:
        values = ' '.join(map(str, array))
        f.write(f"{user_defined_name} {{{values}}}\n")
    
    print(f"Array successfully written to {output_file_path}")
def general_loader(file_path, *columns):
    # Open the file and read data
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('%') or not line.strip() or 'sweep_target' in line:
                continue
            parts = line.split()
            row = []
            for part in parts:
                try:
                    row.append(float(part))
                except ValueError:
                    row.append(part)
            data.append(row)

    # Convert data to numpy array
    data = np.array(data)

    # Extract specified columns
    extracted_columns = []
    for col in columns:
        extracted_columns.append(data[:, col])

    # Return the extracted columns as separate numpy arrays
    return tuple(extracted_columns)
def parse_table(filepath, radius_position, to_table=False):
    def try_convert(value):
        try:
            return int(value)
        except ValueError:
            try:
                return float(value)
            except ValueError:
                return value  
    radius_dict = {}
    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith('%'):
                continue
            data = line.split()
            converted_data = [try_convert(entry) for entry in data]
            try:
                radius = round(float(converted_data[radius_position]))
            except (ValueError, IndexError):
                continue  
            if radius not in radius_dict:
                radius_dict[radius] = []
            radius_dict[radius].append(converted_data)
    if to_table:
        for radius in radius_dict:
            headers = [f"Col{i}" for i in range(len(radius_dict[radius][0]))]
            radius_dict[radius] = pd.DataFrame(radius_dict[radius], columns=headers)
    return radius_dict
def fit_refractive_index_cauchy(wavelength_nm, refractive_index, order=4, 
                                initial_guess=None, wi=300, wf=1600, n_points=800,
                                xi=None, xf=None):
    """
    Fits a refractive index curve using the Cauchy formula within a specified range.
    
    Parameters:
    - wavelength_nm: Array of wavelengths (in nm) for the refractive index data.
    - refractive_index: Array of refractive index values corresponding to wavelength_nm.
    - order: The order of the Cauchy formula to use (2, 4, or 6).
    - initial_guess: Optional initial guess for curve fitting. If None, a default guess is used.
    - wi: Start wavelength (in nm) for the fitted curve output.
    - wf: End wavelength (in nm) for the fitted curve output.
    - n_points: Number of points for the fitted wavelength range.
    - xi: Lower bound of wavelength range (in nm) to be used for fitting. Defaults to min of wavelength_nm.
    - xf: Upper bound of wavelength range (in nm) to be used for fitting. Defaults to max of wavelength_nm.

    Returns:
    - wavelength_fitted_nm: Array of wavelengths from wi to wf, equally spaced with n_points.
    - refractive_index_fitted: The corresponding refractive index values from the fitted Cauchy formula.
    """
    
    # Ensure inputs are NumPy arrays
    wavelength_nm = np.asarray(wavelength_nm)
    refractive_index = np.asarray(refractive_index)
    
    # Validate order
    if order not in [2, 4, 6]:
        raise ValueError("Order must be 2, 4, or 6.")
    
    # Define default initial guess if none is provided
    if initial_guess is None:
        initial_guess = [1.5, 1e-18, 3e-25, 1e-36][:order//2 + 1]

    # Define the Cauchy formula based on order
    def cauchy_formula(lambda_m, A, B, C=0, D=0):
        if order == 2:
            return A + B / lambda_m**2
        elif order == 4:
            return A + B / lambda_m**2 + C / lambda_m**4
        elif order == 6:
            return A + B / lambda_m**2 + C / lambda_m**4 + D / lambda_m**6

    # Apply fitting range if xi and xf are provided
    if xi is None:
        xi = np.min(wavelength_nm)
    if xf is None:
        xf = np.max(wavelength_nm)
    
    # Filter data within the fitting range
    mask = (wavelength_nm >= xi) & (wavelength_nm <= xf)
    wavelength_fit = wavelength_nm[mask]
    refractive_index_fit = refractive_index[mask]

    # Convert wavelength to meters for fitting
    wavelength_fit_m = wavelength_fit * 1e-9

    # Perform the curve fitting
    popt, pcov = curve_fit(cauchy_formula, wavelength_fit_m, refractive_index_fit, p0=initial_guess)

    # Calculate chi-squared
    residuals = refractive_index_fit - cauchy_formula(wavelength_fit_m, *popt)
    chi_squared = np.sum((residuals**2))
    degrees_of_freedom = len(refractive_index_fit) - len(popt)
    reduced_chi_squared = chi_squared / degrees_of_freedom if degrees_of_freedom > 0 else np.nan

    # Print fitted parameters and chi-squared information in scientific notation
    formatted_params = ", ".join([f"{param:.5e}" for param in popt])
    print(f"Fitted Parameters (A, B, C, D): [{formatted_params}]")
    print(f"Reduced Chi-Squared: {reduced_chi_squared:.5e}")

    # Create the specified output wavelength range
    wavelength_fitted_nm = np.linspace(wi, wf, n_points)
    wavelength_fitted_m = wavelength_fitted_nm * 1e-9

    # Calculate the fitted refractive index
    refractive_index_fitted = cauchy_formula(wavelength_fitted_m, *popt)

    return wavelength_fitted_nm, refractive_index_fitted
def quick_bend_Q(file_path,xtick,radius_col=0,a=0.3,TE_TM=None,re_neff_col=4,filter_ratio=True,imag_col=5,ratio_col=6 ,color_opt='blue', mcolor='blue', name='Long', markersize=6,swap_dir=0,tick_font_size=6,xi=0,xf=0,threshold=0.9):
    """
    Plots the bending quality factor (Q) for specified radius values using data from a given file.
    
    Parameters:
    - file_path (str): Path to the data file containing radius, real, and imaginary values.
    - color_opt (str): Color for the main plot lines. Default is 'blue'.
    - mcolor (str): Color for the markers. Default is 'blue'.
    - name (str): Label for the plot legend. Default is 'Long'.
    - size (int): Size of the markers in the plot. Default is 6.
    
    Returns:
    - None. Displays a log-log plot of Bending Q vs. Radius.
    """
    r1 = np.array([])
    q1 = np.array([])
    r2 = np.array([])
    q2 = np.array([])   
    if swap_dir:
        bb = ['$→$','$↑$']
    else:
        bb = ['$↑$', '$→$']
    if TE_TM is None:
        jj = [0, 1]
        bb = bb
    else:
        jj= [TE_TM]
        bb=[bb[TE_TM]]
    # Plot using shared colors for the same Δx
    for rank, marker in zip(jj, bb):
        r, real, imag = get_dog(radius_col, re_neff_col, file_path, filter_ratio=filter_ratio, ratio_col=ratio_col, extra_info_col=imag_col, rank=rank,threshold=threshold)
        q=np.abs(real / imag / 2)
        plt.plot(
            r, q, 
            marker=marker, linestyle='-', label=name, 
            alpha=a, color=color_opt, markersize=markersize, 
            markerfacecolor=mcolor, markeredgecolor=mcolor
        )
        if rank == 0:
            r1 = r
            q1 = q
        else:
            r2=r
            q2=q

    # Add horizontal reference line
    if xi is not 0:
        plt.plot([xi, xf], [10e8, 10e8], linestyle='--', color='black')
    else:
        plt.plot([r1[0], r1[-1]], [10e8, 10e8], linestyle='--', color='black')

    # Configure plot settings
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Radius (µm)')
    plt.ylabel('Bending Q')

    # Set tick font size
    tick_font_size = tick_font_size  # Adjust as needed
    plt.tick_params(axis='both', which='major', labelsize=tick_font_size)
    # Customize x-axis ticks and formatting
    plt.xticks(xtick)
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    formatter.set_useOffset(False)  # Disable offset
    plt.gca().xaxis.set_major_formatter(formatter)
    # Add legend and display plot
    plt.legend(loc='upper left')
    plt.grid(True, which="both", linestyle='--', linewidth=0.5)
    #plt.show()
    return r1,q1,r2,q2
def annotator(r,q,index,label="Q(M)",color='b'):
    for i in index: 
        temp=q[i]/1000000
        if temp>1000:
            info = f'{(temp/1000):.2f}B'
        else:
            info=f'{(q[i]/1000000):.2f}M'
        plt.annotate(info,xy=(r[i],q[i]),xytext=(r[i],q[i]*1.05),ha="right", va="bottom" ,color=color)
def merger(r1, n1, r2, n2):
    """
    Merges two pairs of arrays (r1, n1) and (r2, n2) into a single pair (r, n),
    sorts them by the values in r in ascending order, and maintains the pairings.
    
    Parameters:
    - r1 (array-like): Array of r values in the first group.
    - n1 (array-like): Array of n values paired with r1.
    - r2 (array-like): Array of r values in the second group.
    - n2 (array-like): Array of n values paired with r2.
    
    Returns:
    - r_merged (numpy.ndarray): Merged and sorted array of r values.
    - n_merged (numpy.ndarray): Array of n values corresponding to sorted r_merged.
    """
    
    # Merge the r and n arrays
    r_combined = np.concatenate((r1, r2))
    n_combined = np.concatenate((n1, n2))

    # Sort by r values while keeping the n pairing
    sorted_indices = np.argsort(r_combined)
    r_merged = r_combined[sorted_indices]
    n_merged = n_combined[sorted_indices]

    return r_merged, n_merged
def n(lambda_nm, material_file): 
    data = np.loadtxt(material_file)
    if data[:,0][0] >100:
        wavelength_nm = data[:, 0]
    else:
        wavelength_nm = data[:, 0] * 1e9  

    refractive_index = data[:, 1]
    closest_idx = np.argmin(np.abs(wavelength_nm - lambda_nm))
    return wavelength_nm[closest_idx], refractive_index[closest_idx]
def beta_r(radius_um,offset_um,lambda_nm,core,cladding=1,width_core_um = 0, sweep_r=0):
    ida = lambda_nm*1e-9  
    nCore = n(lambda_nm,core)[1]
    if cladding==1:
        nClad=1
        print('No cladding file received. Assume nClad=1.')
    else:
        nClad = n(lambda_nm,cladding)[1]  
    nCore = n(lambda_nm,core)[1]
    dr = (offset_um+width_core_um/3)*1e-6      
    radius = radius_um*1e-6   

    # Calculate the expression
    inside_sqrt = nClad**2 - (nCore * radius / (radius + dr))**2
    result = ida / cmath.sqrt(inside_sqrt)  # Use cmath.sqrt to handle complex numbers if any
    result_um = result*1e6
    minairoff = (nCore/nClad-1)*radius_um
    print(f"Ida_Radial({lambda_nm}) = {result_um} [um]")
    print(f"minimum air offset is {minairoff:.4f} [um]")

    if sweep_r:
        sweep_radius = np.linspace(0.1,10,100)*radius_um
        def sweep(radius):
            dr = (offset_um+width_core_um/3)*1e-6      
            radius = radius*1e-6   
            inside_sqrt = nClad**2 - (nCore * radius / (radius + dr))**2
            if inside_sqrt < 0:
                result_um = 0
            else:
                result = ida / np.sqrt(inside_sqrt)  # Use np.sqrt for real numbers
                result_um = result * 1e6  # Convert back to micrometers
            return result_um
        plt.figure(figsize=(12,8),dpi=100)
        ida_results = np.array([sweep(r) for r in sweep_radius])
        plt.scatter(sweep_radius,ida_results)
        plt.xlabel('Radius (um)')
        plt.ylabel('Ida_Radial (um)')
def model_loader(lambda_nm,model_in,cladding,core,pml_um = -1, air_offset_um = -1, note=None):

    # set the time stamp
    if note is not None:
        model_in.description('readme',f'{note} @{lambda_nm}[nm]')
    else:
        note = datetime.now().strftime("%m/%d/%Y %I:%M%p")
        try: 
            model_in.description('readme',f'{note} @{lambda_nm}[nm]')
        except ValueError:
            check = model_in.description('readme')
            print(f'readme has been set to: {check}\n')

    # prepare the values 
    if cladding==1:
        nClad = 1
    else: 
        nClad = n(lambda_nm,cladding)[1]
    nCore = n(lambda_nm,core)[1]

    #load them into COMSOL
    model_in.parameter('nCore',f'{nCore}')
    model_in.parameter('nClad',f'{nClad}')
    model_in.parameter('ida',f'{lambda_nm}[nm]')

    if pml_um > 0:
        model_in.parameter('pml_extend',f'{pml_um}[um]')
    if air_offset_um > 0:
        model_in.parameter('air_offset',f'{air_offset_um}[um]')
    #These five below needs to be printed and compare
    ida = model_in.parameter('ida')
    coren = model_in.parameter('nCore',evaluate=1)
    cladn = model_in.parameter('nClad',evaluate=1)
    pml = model_in.parameter('pml_extend')
    air_offset = model_in.parameter('air_offset')
    idaPLM = model_in.parameter('idaPLM',evaluate=1)
    rmodel = model_in.parameter('radius')
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
    print(f'idaPLM({rmodel},{ida}): {idaPLM}')
def export_table_tag(model, table_name, export_path): #This is an intermediate function that export a table
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
def waveguide_resonator(path2d,style='turbo',width=0,height=0,radius=0,qq=1,strong=1, label=None):
    #qq is a parameter that tunes the display
    # Load the data, skipping the metadata lines
    data = np.loadtxt(path2d, skiprows=9)

    # Convert units from meters to micrometers (um)
    R = data[:, 0] * 1e6  # x-axis (R) in micrometers
    Z = data[:, 1] * 1e6  # y-axis (Z) in micrometers
    Color = data[:, 2]  # energy distribution remains unchanged

    # Define grid for interpolation
    R_grid = np.linspace(R.min(), R.max(), 500)
    Z_grid = np.linspace(Z.min(), Z.max(), 500)
    R_grid, Z_grid = np.meshgrid(R_grid, Z_grid)
    Color_grid = griddata((R, Z), Color, (R_grid, Z_grid), method='cubic')
    if 1-strong:
        # Interpolate data onto grid
        info='linear'
        Color_grid=Color_grid
    else:
        # Apply logarithmic scaling to emphasize weaker traces
        Color_grid = np.log(np.abs(Color_grid) + qq)  # small value added to avoid log(0)
        info = 'log'

    a = 1
    # Plot the data with logarithmic scaling
    plt.figure(figsize=(8, 6))
    plt.contourf(R_grid, Z_grid, Color_grid, levels=300, cmap=style)
    plt.colorbar(label=info)

    # Define the hollow square parameters in micrometers
    if width is not 0:

        x_center = radius
        y_center = 0 

        # Create the hollow square with an extended bottom border
        square = patches.Rectangle((x_center - width / 2, y_center - height / 2),
                                width, height, linewidth=a, edgecolor='white', facecolor='none',linestyle='-')
        plt.gca().add_patch(square)

        # Draw the extended line along the bottom edge of the square
        x_left = x_center - width / 2
        x_right =x_center + width / 2
        y_bottom=-height/2
        if label is None:
            plt.plot([R.min(), x_left], [y_bottom, y_bottom], color='white', linewidth=a, linestyle='-',label=f'R={radius}um, {width}x{height}um')
        else: 
            plt.plot([R.min(), x_left], [y_bottom, y_bottom], color='white', linewidth=a, linestyle='-',label=f'R={radius}um, {width}x{height}um, {label}')
        plt.plot([x_right, R.max()], [y_bottom, y_bottom], color='white', linewidth=a, linestyle='-')


    # Label the axes in micrometers
    plt.xlabel("R  [µm]")
    plt.ylabel("Z  [µm]")
    plt.grid(0)
    plt.legend()
    #plt.title(f"TE, E-field Rightward_{style}")
def n_square(wavelength, *args):
    # Convert the wavelength to micrometers if needed
    if wavelength[0] > 100:  # Convert nm to um
        wavelength = wavelength / 1000
    elif wavelength[0] < 0.001:  # Convert m to um
        wavelength = wavelength * 1e6

    # Extract Sellmeier coefficients (assuming 3 terms for B and C)
    b1, c1, b2, c2, b3, c3 = args
    
    # Sellmeier equation
    result = 1 + (b1 * wavelength**2) / (wavelength**2 - c1) \
               + (b2 * wavelength**2) / (wavelength**2 - c2) \
               + (b3 * wavelength**2) / (wavelength**2 - c3)
    return result
def fit_Sellmeier(actual_n, actual_wavelength, b1=0.696166, c1=0.68404e-1, b2=0.4079426, c2=0.1162414, b3=0.89747940, c3=0.989616e1,xi=500,xf=1700,step=1,index_begin=None,index_end=None):
    # Convert the wavelength to micrometers if needed
    if actual_wavelength[0] > 100:  # Convert nm to um
        actual_wavelength = actual_wavelength / 1000
    elif actual_wavelength[0] < 0.001:  # Convert m to um
        actual_wavelength = actual_wavelength * 1e6

    # Initial guess for the fitting parameters
    initial_guess = [b1, c1, b2, c2, b3, c3]
    if index_begin is None:
        # Perform curve fitting to optimize the Sellmeier coefficients
        popt, pcov = curve_fit(n_square, actual_wavelength, actual_n**2, p0=initial_guess,maxfev = 54000)
    else:
        popt, pcov = curve_fit(n_square, actual_wavelength[index_begin:index_end], (actual_n[index_begin:index_end])**2, p0=initial_guess,maxfev = 54000)
    # Generate a range of wavelengths for the fitted curve
    fitted_wavelength = np.arange(xi, xf, step)  # Wavelength range in nm
    # Compute the fitted refractive indices using the optimized coefficients
    fitted_n = np.sqrt(n_square(fitted_wavelength, *popt))

    # Print the fitted Sellmeier coefficients
    print(f"Fitting parameters: b1 = {popt[0]:.6f}, c1 = {popt[1]:.6e}, "
          f"b2 = {popt[2]:.6f}, c2 = {popt[3]:.6e}, "
          f"b3 = {popt[4]:.6f}, c3 = {popt[5]:.6e}")

    return fitted_wavelength, fitted_n,popt
def plotly_plot(x1, y1, x2=None, y2=None, x_label="X-axis", y_label="Y-axis", title="Plot", filename="plot.html"):
    """
    Creates an interactive plot using Plotly with one or two datasets and opens it in the web browser.

    Parameters:
    - x1: array-like, data for the x-axis of the first dataset.
    - y1: array-like, data for the y-axis of the first dataset.
    - x2: array-like, data for the x-axis of the second dataset (optional).
    - y2: array-like, data for the y-axis of the second dataset (optional).
    - x_label: str, label for the x-axis (default is "X-axis").
    - y_label: str, label for the y-axis (default is "Y-axis").
    - title: str, title of the plot (default is "Plot").
    - filename: str, name of the HTML file to save the plot as (default is "plot.html").

    Returns:
    - None (opens the plot in a web browser).
    """
    # Create the plot
    fig = go.Figure()

    # Add first dataset (fitted data or measured data)
    fig.add_trace(go.Scatter(
        x=x1, y=y1,
        mode='markers+lines',
        name='Dataset 1',
        marker=dict(size=5, color='blue'),
        line=dict(width=2)
    ))

    # If a second dataset is provided, add it to the plot
    if x2 is not None and y2 is not None:
        fig.add_trace(go.Scatter(
            x=x2, y=y2,
            mode='markers+lines',
            name='Dataset 2',
            marker=dict(size=5, color='red'),
            line=dict(width=2, dash='dash')  # Dashed line for distinction
        ))

    # Set the axis labels and title
    fig.update_layout(
        title=title,
        xaxis_title=x_label,
        yaxis_title=y_label,
        hovermode="x unified"
    )

    # Add a grid
    fig.update_xaxes(showgrid=True)
    fig.update_yaxes(showgrid=True)

    # Save the plot as an HTML file
    pyo.plot(fig, filename=filename, auto_open=False)

    # Open the saved HTML file in the web browser
    webbrowser.open(filename)
def load_dispersion(file_path, *columns, label_column=None, neff_column=None, filter_column=None, threshold=0.9, rank=1):
    """
    Load data from a file and extract specified columns. Optionally, group data by a label column
    and extract the row with the nth largest neff value within each group, then sort by label_column.

    Additionally, filter rows where values in the filter_column are below the given threshold.

    Parameters:
    - file_path: str, path to the data file.
    - *columns: int, indices of the columns to extract.
    - label_column: int, index of the label column (e.g., 'sweep_target'). Default is None.
    - neff_column: int, index of the neff column. Default is None.
    - filter_column: int, index of the column used for filtering rows. Default is None.
    - threshold: float, the threshold value for filtering rows. Default is 0.9.
    - rank: int, which largest value to return (e.g., 1 for largest, 2 for second largest, etc.)

    Returns:
    - tuple of numpy arrays corresponding to the extracted columns, sorted by label_column.
    """
    # Open the file and read data
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('%') or not line.strip() or 'sweep_target' in line:
                continue
            parts = line.split()
            row = []
            for part in parts:
                try:
                    row.append(float(part))
                except ValueError:
                    row.append(part)  # Keep as string or complex if it cannot be converted to float
            data.append(row)

    # Convert the data to a numpy array (keep as strings/complex numbers where applicable)
    data = np.array(data, dtype=object)  # dtype=object to preserve mixed types

    # Apply filtering based on the filter_column and threshold
    if filter_column is not None:
        data = data[data[:, filter_column].astype(float) >= threshold]

    # If label_column and neff_column are specified, process the data accordingly
    if label_column is not None and neff_column is not None:
        # Extract unique labels (convert only the label_column to float)
        labels = np.unique(data[:, label_column].astype(float))
        # Initialize a list to store rows with the nth largest neff for each label
        max_neff_rows = []
        for label in labels:
            # Get rows corresponding to the current label (again convert only label_column)
            label_rows = data[data[:, label_column].astype(float) == label]
            # Sort the rows based on neff values in descending order
            sorted_rows = label_rows[np.argsort(label_rows[:, neff_column].astype(float))[::-1]]
            # Ensure rank is valid; if it's greater than available rows, pick the last row
            row_index = min(rank - 1, len(sorted_rows) - 1)
            # Select the row corresponding to the nth largest neff
            max_neff_row = sorted_rows[row_index]
            max_neff_rows.append(max_neff_row)
        # Convert the list to a numpy array
        data = np.array(max_neff_rows, dtype=object)

    # Sort the data by the label_column (convert only label_column to float)
    if label_column is not None:
        data = data[np.argsort(data[:, label_column].astype(float))]  # Sort based on the label_column as float

    # Extract specified columns
    extracted_columns = []
    for col in columns:
        extracted_columns.append(np.array(data[:, col], dtype=float))

    # Return the extracted columns as separate numpy arrays
    return tuple(extracted_columns)
def partition_table_by_column(filepath, column_index, output_dir,name):
    """
    Partitions a table based on unique values in a specified column and includes header comments.
    
    Parameters:
    - filepath (str): Path to the input file.
    - column_index (int): Index of the column to partition by (0-based).
    - output_dir (str): Directory where partitioned files will be saved.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Read header lines (lines that start with '%')
    with open(filepath, 'r') as file:
        header_lines = []
        for line in file:
            if line.startswith('%'):
                header_lines.append(line)
            else:
                break  # Stop reading after the header lines

    # Load the table into a DataFrame, skipping the header lines
    data = pd.read_csv(filepath, delim_whitespace=True, skiprows=len(header_lines), header=0)
    
    # Extract the column label from the specified index
    column_label = data.columns[column_index]
    temp=[]
    # Group the data by unique values in the specified column
    for label, group in data.groupby(data.columns[column_index]):
        # Define the output file name based on the unique value in the specified column
        output_file = os.path.join(output_dir, f"{name}_{label}.txt")
        
        # Save the header lines and the group data
        with open(output_file, 'w') as f:
            # Write the header lines
            f.writelines(header_lines)
            
            # Write the partitioned data with tab-separated format, including column headers
            group.to_csv(f, sep='\t', index=False)
        temp.append(label)
        print(f"Saved partition with {name} = {label} to {output_file}")
    print(f'\n\n {temp}')
# Example usage:
# filepath = "path/to/your/input.txt"
# column_index = 0  # Column to partition by, e.g., 'width' in this case is column 0
# output_dir = "path/to/output_directory"
# partition_table_by_column(filepath, column_index, output_dir)
def index_position(raw_headers, column_name):
    """
    Finds the index of the column containing the given substring from a raw string of headers.
    
    Args:
        raw_headers (str): A raw string of column headers.
        column_name (str): The substring to search for in the headers.
    
    Returns:
        int: The index of the column matching the substring, or -1 if not found.
    """
    # Split the raw header string into individual column names
    headers = [header.strip() for header in raw_headers.split("          ") if header.strip()]
    # Iterate to find the matching column
    for index, header in enumerate(headers):
        
        if column_name.lower() in header.lower():
            return index
    return -1  # Return -1 if the column is not found