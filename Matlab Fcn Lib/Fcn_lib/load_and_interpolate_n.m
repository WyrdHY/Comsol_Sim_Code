function n_interp = load_and_interpolate_n(filepath)
% Load the data from the provided file path
    data = load(filepath);
    wavelength_nm = data(:,1);  % First column: wavelength (nm)
    n_values = data(:,2);       % Second column: refractive index (n)

    % Create the interpolation function
    n_interp = @(wavelength) interp1(wavelength_nm, n_values, wavelength, 'spline');

    % Display a confirmation message
    disp('Interpolation function created successfully.');
    disp(['Wavelength range: ', num2str(min(wavelength_nm)), ' nm to ', num2str(max(wavelength_nm)), ' nm']);
end

%% Gaussian Mode
