%% Initialization
addpath('C:\Program Files\COMSOL\COMSOL55\Multiphysics\mli')
import com.comsol.model.*
import com.comsol.model.util.*
format long;
mphstart(2036); 
%%
modelPath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Lensed_Fiber_to_Waveguide_Coupling\Straight_Wg.mph";
savePath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Lensed_Fiber_to_Waveguide_Coupling\Straight_Wg.mph";

model = mphload(modelPath);   
ModelUtil.showProgress(true);

%% Navigator
mphnavigator(model);

%% Configure Parameters for Dispersion Single Ring
fpath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Material\Sellmeier Fitting\Sellmeier_2%_Long Anneal.txt";
n = load_and_interpolate_n(fpath);
core_width = 12;
core_height = 2;
ida = 532;
n_initial_guess = n(ida);
disp(n_initial_guess)
%%
model.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
model.param.set('initial_guess', num2str(n_initial_guess, '%.5f'));
model.param.set('core_width', [num2str(core_width, '%.4f'),'[um]']);
model.param.set('core_height', [num2str(core_height, '%.4f'),'[um]']);
disp('Configured');
%% Run Study For TE mode
indicator = flag(333); %Do ewfd
studyName = indicator.std; 
disp('Running')

model.study(studyName).run();
%%
neff = mphglobal(model, 'real(ewfd.neff)', 'dataset', indicator.dset);
[~, eigs] = sort(neff, 'descend');
top2Indices = eigs(1:2);

temp1 = mphglobal(model, '(ez)', 'dataset', indicator.dset,'solnum',top2Indices(1));
temp2 = mphglobal(model, '(ez)', 'dataset', indicator.dset,'solnum',top2Indices(2));

if temp1>temp2
    iTE = top2Indices(1);
    iTM = top2Indices(2);
else
    iTE = top2Indices(2);
    iTM = top2Indices(1);
end
disp('Body Momentum Done')

modes = {'TE', 'TM'};
solnums = [iTE, iTM];  
titles = {'TE Mode, E-Up', 'TM Mode, E-Right'};
figure;
for i = 1:2
    model.result('pg1').set('data', indicator.dset);
    model.result('pg1').set('solnum', num2str(solnums(i)));
    model.result('pg1').feature('surf1').set('expr', "ewfd.normE");  
    model.result('pg1').feature('surf1').set('descr', '|E|');        
    subplot(1, 2, i);
    mphplot(model, 'pg1');
    title(titles{i});
end
drawnow;
disp(['Study ', studyName, ' completed successfully.']);
%% Loop through different beam waist
% Prepare Mesh Grid
w0_list = (0.5:0.2:5)*1e-6;
polar_list = [0,1];
eta_list = zeros(2,length(w0_list));

for polar = polar_list
    i=1;
    for w0 = w0_list
        % 0-iTM, 1-iTE
        if polar ==1 
            solnum=iTE;
            j = 1;
        else
            solnum=iTM;
            j=2;
        end
        %Para for Gaussian Beam
        nx = 500;  
        ny = 300;  
        %Meshing 
        rx = max(2*w0,10e-6);
        ry = max(2*w0,10e-6);
        xVec = linspace(-rx, rx, nx);
        yVec = linspace(-ry, ry, ny);
        [X, Y] = meshgrid(xVec, yVec); 
        coords = [X(:)'; Y(:)'];
        %Extract
        Ex_data = reshape(mphinterp(model, 'ewfd.Ex','coord',coords,'dataset',indicator.dset,'solnum',solnum), ny, nx);
        Ey_data = reshape(mphinterp(model, 'ewfd.Ey','coord',coords,'dataset',indicator.dset,'solnum',solnum), ny, nx);
        %Generate Gaussian 
        [E_gauss_x, E_gauss_y] = Gauss(w0,X,Y,polar); 
        % Numerically Calculate the Integral
        dx = xVec(2)-xVec(1);
        dy = yVec(2)-yVec(1);
        I12_mat = conj(Ex_data).*E_gauss_x + conj(Ey_data).*E_gauss_y;
        I12 = sum(I12_mat(:))*dx*dy;
        I11 = sum((abs(Ex_data).^2 + abs(Ey_data).^2), 'all')*dx*dy;
        I22 = sum((abs(E_gauss_x).^2 + abs(E_gauss_y).^2), 'all')*dx*dy;
        eta_list(j,i) = abs(I12)/sqrt(I11*I22);
        i=i+1;
    end
end
%% Plot the beam_waist stats
w0_list_um = w0_list * 1e6; 
te_color = [0 0.447 0.741];   % Nice blue color for TE
tm_color = [0.850 0.325 0.098]; % Nice red-orange color for TM

figure;
hold on;

plot(w0_list_um, eta_list(1,:), '--d', ...
    'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', te_color, 'Color', te_color, 'DisplayName', 'TE Mode');

plot(w0_list_um, eta_list(2,:), '--^', ...
    'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', tm_color, 'Color', tm_color, 'DisplayName', 'TM Mode');

xlabel('Beam Waist w_0 [\mum]', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Coupling Efficiency \eta', 'FontSize', 16, 'FontWeight', 'bold');
title('Coupling Efficiency vs Beam Waist', 'FontSize', 18, 'FontWeight', 'bold');

grid on;
set(gca, 'GridLineStyle', '--', 'FontSize', 14);

legend('Location', 'northwest', 'FontSize', 14);

hold off;

%% Check the Plot
compare_plot(2,iTE,1,model,indicator);

%% Saving
mphsave(model, savePath);
disp(['Model saved to: ', savePath]);
%% Exit 
% ModelUtil.disconnect();
ModelUtil.remove('model');
disp('Model removed from memory.');

%%
function output = flag(x)
   % x = 0, means small region, ewfd2
   % x = 333 means no need to distinguish them
   output = struct();
   if x==0
       output.dset = 'dset1';
       output.normE = 'ewfd2.normE';
       output.std = 'std1';
   else
       output.dset = 'dset2';
       output.normE = 'ewfd1.normE';      
       output.std = 'std2';
   end 

   if x == 333
        output.dset = 'dset1';
        output.normE = 'ewfd1.normE';      
        output.std = 'std1';
   end
end

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

function [E_gauss_x,E_gauss_y]=Gauss(w0,X,Y,polarization)
    % polarization = 1 means TE, Eup
    % polarization = 0 means TM, Eright
    x0 = 0; y0 = 0; A  = 1;  
    if polarization
        E_gauss_y = A * exp(-((X - x0).^2 + (Y - y0).^2)/(w0^2));
        E_gauss_x = zeros(size(E_gauss_y));
    else
        E_gauss_x = A * exp(-((X - x0).^2 + (Y - y0).^2)/(w0^2));
        E_gauss_y = zeros(size(E_gauss_x));
    end
end

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