%% Include the Function Lib
basefolder = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Matlab Fcn Lib\Fcn_lib";
addpath(genpath(basefolder)); 
%% Initialization
addpath('C:\Program Files\COMSOL\COMSOL55\Multiphysics\mli')
import com.comsol.model.*
import com.comsol.model.util.*
format long;
mphstart(2036); 
%%
modelPath1 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Pulley_Coupler_2025_05_16\Single Ring.mph";
savePath1 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Pulley_Coupler_2025_05_16\Single Ring_Save.mph";


model1 = ModelUtil.create('Model1');
model1 = mphload(modelPath1, 'Model1');

ModelUtil.showProgress(true);
%% Navigator
mphnavigator(model1);
    % Std1 is for small Region ewfd2
    % Std2 is for large Region ewfd
%% Global Parameters
ida=1064;
radius=3;

n=load_and_interpolate_n("C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\measured_refractive_index\Ge-SiO2 (2% doping RTA) (Shijia-UCSB 2024a n 0.459-1.688 µm).txt");
nCore = n(ida);
nClad = 1;
betaR(nClad, nCore, ida, 5, air_offset, 0, 1)


%%
% Configure Parameters for Body Ring
%{
model1.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
model1.param.set('pml_extend', [num2str(pml_extend, '%.2f'),'[um]']);
model1.param.set('radius', [num2str(radius, '%.3f'),'[um]']);
model1.param.set('r_', [num2str(r_, '%.3f'),'[um]']);
model1.param.set('initial_guess', ['nCore*(radius+r_)']);
model1.param.set('air_offset', [num2str(air_offset, '%.4f'),'[um]']);
disp('Body Configured')
%}

%% D1 D2 for the body ring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ida = 1064;
width = 12;
height = 2;
radius = 1000;
core_mesh_size= 100;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocate Array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = 5;
wavelength_range = ida-A:1:ida+A;
num_points = length(wavelength_range);
wavelength = zeros(1, num_points); % Wavelengths in meters
omega = zeros(1, num_points);  % Omega values
mu = zeros(1, num_points);        % Beta values
j=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startTime = datetime('now');  % Record start time
disp(['Start time: ', datestr(startTime, 'yyyy-mm-dd HH:MM:SS')]);
for x = wavelength_range
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Configure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model1.param.set('ida', [num2str(x, '%.3f'),'[nm]']);
    model1.param.set('radius', [num2str(radius, '%.3f'),'[um]']);
    model1.param.set('height', [num2str(height, '%.3f'),'[um]']);
    model1.param.set('width', [num2str(width, '%.3f'),'[um]']);
    model1.param.set('initial_guess', ['nCore*(radius+0.85*0.5*width)']);
    model1.param.set('core_mesh_size', [num2str(core_mesh_size, '%.3f'),'[nm]']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run Simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    [iTE,iTM,indicator,skip] = calculate(model1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hFig=figure;
    log_plot = 1;
    arrow = 1;
    tinfo = sprintf("j=%d;log plot=%d",j,log_plot);
    visualize(model1,indicator,iTE,log_plot,arrow,tinfo)
    % Save Temporal Photo In Case You Need It Later
    outFolder = 'C:\Users\Dirk\Desktop\Temp_Output';
    filename = fullfile(outFolder, tinfo+ '.fig');
    savefig(hFig, filename);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract Dispersion Information
        % w(mu) = D0 + D1* + D2* ....
        % mu is beta_eff=beta*R
        % omega is just omega
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ewfd_omega = mphglobal(model1, '(ewfd.omega)', 'dataset', indicator.dset,'solnum',iTE);
    ewfd_beta_eff = mphglobal(model1, '(ewfd.beta)', 'dataset', indicator.dset,'solnum',iTE);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    omega(j)=ewfd_omega;
    mu(j) = ewfd_beta_eff;
    disp(sprintf('%d/%d',j,num_points));
    j=j+1;
end
endTime = datetime('now');  % Record end time
disp(['End time: ', datestr(endTime, 'yyyy-mm-dd HH:MM:SS')]);
elapsedTime = endTime - startTime;
disp(['Elapsed time: ', char(elapsedTime)]);
%% Fit the Dispersion
c = 299792458;
omega0 = 2*pi*c/ida*1e9;

mu0 = interp1(omega, mu, omega0);

%Taylor expand w(mu) around w0,u0
w = omega-omega0;
u = mu-mu0;

n = 5;
p = polyfit(u, w, n);
u_fit = linspace(min(u), max(u), 300);    
w_fit = p(end) * ones(size(u_fit));
for i = 1:n
    ai = p(end - i);
    w_fit = w_fit + ai * u_fit.^i;
end

D1 = p(end-1);
D2 = 2*p(end-2); %Negative = Normal;Positive = Anomalous

info = sprintf('%d-th order:\nFSR=D1/2pi=%.2fGhz\nD2=%.2f',n,D1/10^9/2/pi,D2);
disp(info);

give_plot = 1;
if give_plot
figure; hold on;
plot(u, w-D1*u, 'o', 'MarkerSize',6, 'DisplayName','data');         
plot(u_fit, w_fit-D1*u_fit, '-', 'LineWidth',2 );
hold off;

xlabel('\Delta\mu');
ylabel('\Delta\omega');
title(info);
grid on;
end
%% Momentum Mismatch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ida = 1064;
width = 12;
height = 2;
radius = 1000;

radius2 = 1008.10495818933;
width2 = 3.00991637866273;

core_mesh_size= 50;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocate Array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1 = 0;
beta2 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startTime = datetime('now');  % Record start time
disp(['Start time: ', datestr(startTime, 'yyyy-mm-dd HH:MM:SS')]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure Body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model1.param.set('ida', [num2str(x, '%.3f'),'[nm]']);
model1.param.set('radius', [num2str(radius, '%.3f'),'[um]']);
model1.param.set('height', [num2str(height, '%.3f'),'[um]']);
model1.param.set('width', [num2str(width, '%.3f'),'[um]']);
model1.param.set('initial_guess', ['nCore*(radius+0.85*0.5*width)']);
model1.param.set('core_mesh_size', [num2str(core_mesh_size, '%.3f'),'[nm]']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
[iTE,iTM,indicator,skip] = calculate(model1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFig=figure;
log_plot = 1;
arrow = 1;
tinfo = sprintf("Body;log plot=%d",log_plot);
visualize(model1,indicator,iTE,log_plot,arrow,tinfo)
outFolder = 'C:\Users\Dirk\Desktop\Temp_Output';
filename = fullfile(outFolder, tinfo+ '.fig');
savefig(hFig, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b1 = mphglobal(model1, '(ewfd.beta)', 'dataset', indicator.dset,'solnum',iTE);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure Bus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model1.param.set('ida', [num2str(x, '%.3f'),'[nm]']);
model1.param.set('radius', [num2str(radius2, '%.3f'),'[um]']);
model1.param.set('height', [num2str(height, '%.3f'),'[um]']);
model1.param.set('width', [num2str(width2, '%.3f'),'[um]']);
model1.param.set('initial_guess', ['nCore*(radius+0.85*0.5*width)']);
model1.param.set('core_mesh_size', [num2str(core_mesh_size, '%.3f'),'[nm]']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
[iTE,iTM,indicator,skip] = calculate(model1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFig=figure;
log_plot = 1;
arrow = 1;
tinfo = sprintf("Bus;log plot=%d",log_plot);
visualize(model1,indicator,iTE,log_plot,arrow,tinfo)
outFolder = 'C:\Users\Dirk\Desktop\Temp_Output';
filename = fullfile(outFolder, tinfo+ '.fig');
savefig(hFig, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2 = mphglobal(model1, '(ewfd.beta)', 'dataset', indicator.dset,'solnum',iTE);



mismatch = b1-b2;


endTime = datetime('now');  % Record end time
disp(['End time: ', datestr(endTime, 'yyyy-mm-dd HH:MM:SS')]);
elapsedTime = endTime - startTime;
disp(['Elapsed time: ', char(elapsedTime)]);















%% Visualize a mode(Debug)
figure;
model1.result('pg1').set('data', indicator.dset);
model1.result('pg1').set('solnum', num2str(iTE));
model1.result('pg1').feature('surf1').set('expr', '(ewfd.normE)'); % Ensure 'surf1' exists in pg1
model1.result('pg1').feature('surf1').set('descr', '|E|');  
model1.result('pg1').feature('arws1').set('expr', {'abs(ewfd.Er)','abs(ewfd.Ez)'});
mphplot(model1, 'pg1');
drawnow;

%% Saving
if 1
    mphsave(model1, savePath1);
    disp(['Model saved to: ', savePath1]);
    ModelUtil.remove('model1');
end
%%
function betaR(n_Clad, n_Core, lambda_nm, radius_um, offset_um, width_core_um, sweep)
    % betaR calculates the radial parameter (ida_Radial) and, if requested,
    % sweeps over radius values to plot its behavior.
    %
    % Parameters:
    %   n_Clad       - Cladding refractive index. (If no cladding is provided, use 1.)
    %   n_Core       - Core refractive index.
    %   lambda_nm    - Wavelength in nanometers.
    %   radius_um    - Fiber radius in micrometers.
    %   offset_um    - Air offset in micrometers.
    %   width_core_um- (Optional) Core width in micrometers (default = 0).
    %   sweep        - (Optional) Boolean flag. If true, perform a sweep and plot (default = false).
    %
    % Example:
    %   betaR(1, 1.45, 1550, 50, 5, 0, true)
    
    % Set default values for optional parameters
    if nargin < 7
        sweep = false;
    end
    if nargin < 6
        width_core_um = 0;
    end

    % Convert lambda from nanometers to meters
    ida = lambda_nm * 1e-9;
    
    % Convert radius from micrometers to meters and compute additional offset dr
    radius_m = radius_um * 1e-6;
    dr = (offset_um + width_core_um/3) * 1e-6;
    
    % Compute the term inside the square root:
    inside_sqrt = n_Clad^2 - ( n_Core * radius_m / (radius_m + dr) )^2;
    
    % Calculate ida_Radial.
    % Note: sqrt of a negative number returns a complex result in Matlab.
    result = ida / sqrt(inside_sqrt);
    % Convert the result back to micrometers
    result_um = result * 1e6;
    % Compute the minimum air offset
    minairoff = (n_Core/n_Clad - 1) * radius_um;
    % Print the results
    fprintf('Ida_Radial(%f nm) = %g [um]\n', lambda_nm, result_um);
    fprintf('Minimum air offset is %g [um]\n', minairoff);
    if sweep
        sweep_radius = linspace(0.1, 10, 100) * radius_um;
        ida_results = zeros(size(sweep_radius));
        for i = 1:length(sweep_radius)
            r = sweep_radius(i);
            r_m = r * 1e-6;
            inside_sqrt_sweep = n_Clad^2 - ( n_Core * r_m / (r_m + dr) )^2;
            if inside_sqrt_sweep < 0
                ida_results(i) = 0;
            else
                ida_results(i) = ida / sqrt(inside_sqrt_sweep) * 1e6;
            end
        end
        figure;
        scatter(sweep_radius, ida_results, 'filled');
        xlabel('Radius (um)');
        ylabel('Ida Radial (um)');
    end
end

function [iTE,iTM,skip]=eigen_locate(model1,indicator)
    skip = 0;
    neff = mphglobal(model1, 'real(ewfd2.neff)', 'dataset', indicator.dset);
    core_total = mphglobal(model1, 'eigen_core_energy/eigen_total_energy', 'dataset', indicator.dset);
    % Filter out modes with core_total < 0.6
    validIndices = find(core_total >= 0.4);
    % Sort the valid eigenvalues in descending order
    [~, sortOrder] = sort(neff(validIndices), 'descend');
    try
    top2Indices = validIndices(sortOrder(1:2));
    temp1 = mphglobal(model1, '(ez2)', 'dataset', indicator.dset,'solnum',top2Indices(1));
    temp2 = mphglobal(model1, '(ez2)', 'dataset', indicator.dset,'solnum',top2Indices(2));
    if temp1>temp2
        iTM = top2Indices(1);
        iTE = top2Indices(2);
    else
        iTM = top2Indices(2);
        iTE = top2Indices(1);
    end
    catch
        iTE=0;
        iTM=0;
        skip = 1; 
    end
end

function [iTE,iTM,indicator,skip]=calculate(model1)
    %Make sure to include the file defining indicator
    indicator = flag(0); % Small Region ewfd2
    studyName = indicator.std; 
    model1.study('std1').feature('mode').set('neigs', '4');
    model1.study(studyName).run();
    [iTE,iTM,skip] = eigen_locate(model1,indicator);

    if skip==1
        iTE=0;
        iTM=0;
        disp("Skipped")
        return 
    end

    neff = mphglobal(model1, 'real(ewfd2.neff)', 'dataset', indicator.dset);
    initial_guess = neff(iTE);
    
    
    model1.param.set('initial_guess', num2str(initial_guess));
    indicator = flag(1); % Large Region ewfd1
    studyName = indicator.std; 
    model1.study('std2').feature('mode').set('neigs', '4');
    model1.study(studyName).run();

    neff = mphglobal(model1, 'real(ewfd.neff)', 'dataset', indicator.dset);
    core_total = mphglobal(model1, 'core_energy/total_energy', 'dataset', indicator.dset);
    validIndices = find(core_total >= 0.2);
    watchdog = 1;

    if numel(validIndices) < 2
        disp(['Skipping ' ' Careful with this']);
        iTE = validIndices(1);
        iTM = 0;
        xxx = mphglobal(model1, '(ex)', 'dataset', indicator.dset,'solnum',iTE);
        watchdog = 0;
    end
    
    if watchdog
        % Sort the valid eigenvalues in descending order
        [~, sortOrder] = sort(neff(validIndices), 'descend');
        top2Indices = validIndices(sortOrder(1:2));
    
        temp1 = mphglobal(model1, '(ez)', 'dataset', indicator.dset,'solnum',top2Indices(1));
        temp2 = mphglobal(model1, '(ez)', 'dataset', indicator.dset,'solnum',top2Indices(2));
        if temp1>temp2
            iTM = top2Indices(1);
            iTE = top2Indices(2);
        else
            iTM = top2Indices(2);
            iTE = top2Indices(1);
        end
    end
end
function [iTE1,iTE2,indicator,skip]=calculate2(model2) %This is for DUAL RING
    model1=model2;
    indicator = flag(0); % Small Region ewfd2
    studyName = indicator.std; 
    model1.study('std1').feature('mode').set('neigs', '4');
    model1.study(studyName).run();
    [iTE,iTM,skip] = eigen_locate(model1,indicator);

    if skip==1
        iTE=0;
        iTM=0;
        disp("Skipped")
        return 
    end

    neff = mphglobal(model1, 'real(ewfd2.neff)', 'dataset', indicator.dset);
    initial_guess = neff(iTE);
    
    
    model1.param.set('initial_guess', num2str(initial_guess));
    indicator = flag(1); % Large Region ewfd1
    studyName = indicator.std; 
    model1.study('std2').feature('mode').set('neigs', '4');
    model1.study(studyName).run();

    neff = mphglobal(model1, 'real(ewfd.neff)', 'dataset', indicator.dset);
    core_total = mphglobal(model1, 'core_energy/total_energy', 'dataset', indicator.dset);
    ex = mphglobal(model1, 'ex', 'dataset', indicator.dset);
    ez = mphglobal(model1, 'ez', 'dataset', indicator.dset);

    % ——— 2) Filter 1: core_total ≥ 0.2 ———
    valid1  = find(core_total >= 0.2);
    
    % ——— 3) Filter 2: ex/ez ≥ 1 ———
    ratios = ex(valid1) ./ ez(valid1);
    valid2 = valid1(ratios >= 1);


    % ——— 4) Select top-2 by ex descending ———
    [~, sortOrder] = sort(ez(valid2), 'descend');
    top2 = valid2(sortOrder(1:min(2,numel(sortOrder))));

    % Assign your mode indices:
    if numel(top2) >= 2
        iTE1 = top2(1);
        iTE2 = top2(2);
    elseif numel(top2) == 1
        iTE1 = top2(1);
        iTE2 = 0;    % only one valid mode
    else
        iTE1 = 0;
        iTE2 = 0;
    end
end
function visualize(model1,indicator,solnum,log_plot,arrow,tinfo)
    
    
    model1.result('pg1').set('data', indicator.dset);
    model1.result('pg1').set('solnum', solnum);
    model1.result('pg1').set('edgecolor', 'white');
    model1.result('pg1').feature('surf1').set('colortable', 'Rainbow');

 
    if log_plot
    model1.result('pg1').feature('surf1').set('expr', "log(ewfd.normE)");  
    model1.result('pg1').feature('surf1').set('rangecoloractive', 'on');
    model1.result('pg1').feature('surf1').set('rangecolormin', '-4');
    model1.result('pg1').feature('surf1').set('rangecolormax', '5.5');
    model1.result('pg1').feature('surf1').set('descr', '|E|');
    else
    model1.result('pg1').feature('surf1').set('expr', "(ewfd.normE)");  
    model1.result('pg1').feature('surf1').set('rangecoloractive', 'off');
    end

    model1.result('pg1').feature('arws1').active(arrow);
    model1.result('pg1').feature('arws1').set('expr', {'abs(ewfd.Er)','abs(ewfd.Ez)'});
    TE_field = mphplot(model1, 'pg1','rangenum',1);

    title(tinfo, 'Interpreter', 'none')
    drawnow;


end