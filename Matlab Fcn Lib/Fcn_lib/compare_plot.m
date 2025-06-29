function compare_plot(w0, solnum, Gauss_polar, model, indicator)
% Function to compare Gaussian beam with COMSOL field
    % w0: Beam waist in um
    % solnum: Solution number (TE or TM mode)
    % Gauss_polar: Polarization direction (0 for x, 1 for y)
    % model: COMSOL model object
    % indicator: Dataset structure
    
    % Parameters for Gaussian Beam
    nx = 500;  
    ny = 300;  
    w0 = w0*1e-6;
    % Meshing
    rx = max(2*w0, 10e-6);
    ry = max(2*w0, 10e-6);
    xVec = linspace(-rx, rx, nx);
    yVec = linspace(-ry, ry, ny);
    [X, Y] = meshgrid(xVec, yVec); 
    coords = [X(:)'; Y(:)'];
    
    % Extract COMSOL fields
    Ex_data = reshape(mphinterp(model, 'ewfd.Ex', 'coord', coords, 'dataset', indicator.dset, 'solnum', solnum), ny, nx);
    Ey_data = reshape(mphinterp(model, 'ewfd.Ey', 'coord', coords, 'dataset', indicator.dset, 'solnum', solnum), ny, nx);

    % Generate Gaussian beam field
    [E_gauss_x, E_gauss_y] = Gauss(w0, X, Y, Gauss_polar); 

    % Convert x and y to micrometers
    xVec_um = xVec * 1e6;
    yVec_um = yVec * 1e6;

    % Create figure with side-by-side subplots
    figure;
    colormap jet;  

    % Plot Gaussian field
    % Set dynamic title based on polarization
    if Gauss_polar == 1
        title_str = sprintf('w0 = %.2f \x03BCm, TE Mode', w0 * 1e6);
    else 
        title_str = sprintf('w0 = %.2f \x03BCm, TM Mode', w0 * 1e6);
    end
    subplot(1,2,1);
    imagesc(xVec_um, yVec_um, sqrt(E_gauss_y.^2+E_gauss_x.^2));
    
    title(title_str);
    set(gca, 'YDir', 'normal');  
    axis image;                 
    colorbar;
    xlabel('x [\mum]', 'FontSize', 12);
    ylabel('y [\mum]', 'FontSize', 12);

    % Plot COMSOL field
    subplot(1,2,2);
    if Gauss_polar == 1
        imagesc(xVec_um, yVec_um, abs(Ey_data));
    else 
        imagesc(xVec_um, yVec_um, abs(Ex_data));
    end
    set(gca, 'YDir', 'normal');  
    title('COMSOL Simulation')
    axis image;                  
    colorbar;




    xlabel('x [\mum]', 'FontSize', 12);
    ylabel('y [\mum]', 'FontSize', 12);

    % Overall figure title
    sgtitle("Compare")

end
%% Dispersion
