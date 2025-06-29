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
modelPath1 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Bending_Q\Bending_Q.mph";
savePath1 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Bending_Q\Bending_Q__2.mph";


model1 = ModelUtil.create('Model1');
model1 = mphload(modelPath1, 'Model1');

ModelUtil.showProgress(true);
%% Navigator
mphnavigator(model1);
    % Std1 is for small Region ewfd2
    % Std2 is for large Region ewfd
%% Global Parameters
ida=1064;
radius=500;
air_offset = 20;
r=80;
width = 10;
pml_extend = 3;
height = 2; 
meshsize = 100;
air_left = 15;

%n=load_and_interpolate_n("C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\Sellmeier Fitting\Sellmeier_4%_PX332_FHD.txt");
n=load_and_interpolate_n("C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\measured_refractive_index\Ge-SiO2 (2% doping RTA) (Shijia-UCSB 2024a n 0.459-1.688 µm).txt");
nCore = n(ida);

%n = load_and_interpolate_n("C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\measured_refractive_index\ICP-PECVD Long Anneal.txt");
n = load_and_interpolate_n("C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\measured_refractive_index\ReMeasured_Good_THOX.txt");

nClad = n(ida);
betaR(nClad, nCore, ida, r, air_offset, width, 1)


%%
% Configure Parameters for Body Ring
model1.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
model1.param.set('nClad', num2str(nClad, '%.5f'));
model1.param.set('nCore', num2str(nCore, '%.5f'));
model1.param.set('pml_extend', [num2str(pml_extend, '%.2f'),'[um]']);
model1.param.set('radius', [num2str(radius, '%.3f'),'[um]']);
model1.param.set('width', [num2str(width, '%.4f'),'[um]']);
model1.param.set('height', [num2str(height, '%.4f'),'[um]']);
model1.param.set('initial_guess', ['nCore*(radius+0[um])']);
model1.param.set('air_offset', [num2str(air_offset, '%.4f'),'[um]']);
model1.param.set('air_left', [num2str(air_left, '%.4f'),'[um]']);
model1.param.set('core_mesh_size', [num2str(meshsize),'[nm]']);
disp('Body Configured')
%% Calculate Bending Q

radius_list = 100:10:200;
result = zeros(5,length(radius_list)); %r,TE,Q TE, TM,Q TM
result(1,:) = radius_list';
i=1;
for r = radius_list
    model1.param.set('radius', [num2str(r, '%.3f'),'[um]']);
    % Prepare Macro for Simulation 
    model1.param.set('initial_guess', ['nCore*(radius+width/2)']);
    indicator = flag(0); % Small Region ewfd2
    studyName = indicator.std; 
    model1.study('std1').feature('mode').set('neigs', '8');
    disp(sprintf('i = %d, radius = %d',i,r));
    model1.study(studyName).run();
    neff = mphglobal(model1, 'real(ewfd2.neff)', 'dataset', indicator.dset);
    [~, eigs] = sort(neff, 'descend');
    top2Indices = eigs(1:2);
    temp1 = mphglobal(model1, '(ez2)', 'dataset', indicator.dset,'solnum',top2Indices(1));
    temp2 = mphglobal(model1, '(ez2)', 'dataset', indicator.dset,'solnum',top2Indices(2));
    
    if temp1>temp2
        a = top2Indices(1);
        b = top2Indices(2);
    else
        a = top2Indices(2);
        b = top2Indices(1);
    end
    initial_guess = neff(b);
    
    
    % Prepare Macro for detailed Scan
    model1.param.set('initial_guess', num2str(initial_guess));
    indicator = flag(1); % Large Region ewfd1
    studyName = indicator.std; 
    model1.study('std2').feature('mode').set('neigs', '8');
    model1.study(studyName).run();
    neff = mphglobal(model1, 'real(ewfd.neff)', 'dataset', indicator.dset);

    [~, eigs] = sort(neff, 'descend');
    top2Indices = eigs(1:2);
    temp1 = mphglobal(model1, '(ez)', 'dataset', indicator.dset,'solnum',top2Indices(1));
    temp2 = mphglobal(model1, '(ez)', 'dataset', indicator.dset,'solnum',top2Indices(2));
    
    if temp1>temp2
        iTM = top2Indices(1);
        iTE = top2Indices(2);
    else
        iTM = top2Indices(2);
        iTE = top2Indices(1);
    end

    check = mphglobal(model1, 'core_energy/total_energy', 'dataset', indicator.dset);

    if check(iTE)>0.3
        beta_TE = mphglobal(model1, 'real(ewfd.beta)', 'dataset', indicator.dset,'solnum',iTE);
        re = mphglobal(model1,'real(ewfd.neff)', 'dataset', indicator.dset,'solnum',iTE);
        im = mphglobal(model1,'imag(ewfd.neff)', 'dataset', indicator.dset,'solnum',iTE);
        Q_TE = re/(im*2);
    else
        beta_TE = 0;
        Q_TE = 0;
    end

    if check(iTM)>0.3
        beta_TM = mphglobal(model1, 'real(ewfd.beta)', 'dataset', indicator.dset,'solnum',iTM);
        re = mphglobal(model1,'real(ewfd.neff)', 'dataset', indicator.dset,'solnum',iTM);
        im = mphglobal(model1,'imag(ewfd.neff)', 'dataset', indicator.dset,'solnum',iTM);
        Q_TM = re/(im*2);
    else
        beta_TM = 0;
        Q_TM = 0;
    end
    result(2,i) = beta_TE;
    result(3,i) = Q_TE;
    result(4,i) = beta_TM;
    result(5,i) = Q_TM;
    disp(i/length(radius_list));
    i = i+1;
end
disp('Done')
%%
im = mphglobal(model1, 'imag(ewfd.beta)', 'dataset', indicator.dset, 'solnum', iTE);
disp(im)

%%
% Extract variables for clarity
r      = result(1,:);   % r values
betaTE = result(2,:);
Q_TE   = abs(result(3,:));
betaTM = result(4,:);
Q_TM   = abs(result(5,:));

figure;

subplot(2,1,1);
% Plot Q_TE vs r using semilogy for log-scale on y-axis
semilogy(r, Q_TE, '--d', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
hold on;
% Plot Q_TM vs r
semilogy(r, Q_TM, '--d', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
% Draw horizontal dashed line at y = 1e9 (no label)
yline(1e9, '--', 'LineWidth', 1);
hold off;
xlabel('r');
ylabel('Q');
title('r vs Q (Log Scale)');
legend('TE','TM','Location','best');
subplot(2,1,2);
plot(r, betaTE, '--d', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
hold on;
plot(r, betaTM, '--d', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
hold off;

xlabel('r');
ylabel('\beta');
title('r vs \beta');
legend('TE','TM','Location','best');

%% Visualize it
figure;
model1.result('pg1').set('data', indicator.dset);
model1.result('pg1').set('solnum', num2str(iTE));
model1.result('pg1').feature('surf1').set('expr', '(ewfd.normE)'); % Ensure 'surf1' exists in pg1
model1.result('pg1').feature('surf1').set('descr', '|E|');  
model1.result('pg1').feature('arws1').set('expr', {'abs(ewfd.Er)','abs(ewfd.Ez)'});
mphplot(model1, 'pg1');
drawnow;
%%
result = zeros(10,5);
disp(result(1,:));

%% Saving
if 1

    mphsave(model1, savePath1);
    disp(['Model saved to: ', savePath1]);
    % Exit 
    % ModelUtil.disconnect();
    ModelUtil.remove('model1');

    disp('Model removed from memory.');
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


%{

% Prepare Macro for Simulation 
model1.param.set('initial_guess', ['nCore*(radius+0[um])']);
indicator = flag(0); % Small Region ewfd2
studyName = indicator.std; 
model1.study('std1').feature('mode').set('neigs', '2');
disp('Running Body Momentum')
model1.study(studyName).run();
neff = mphglobal(model1, 'real(ewfd2.neff)', 'dataset', indicator.dset);
[initial_guess,iMax] = max(neff);


% Prepare Macro for detailed Scan
model1.param.set('initial_guess', num2str(initial_guess));
indicator = flag(1); % Large Region ewfd1
studyName = indicator.std; 
model1.study('std2').feature('mode').set('neigs', '2');
model1.study(studyName).run();
[~, eigs] = sort(neff, 'descend');
top2Indices = eigs(1:2);
temp1 = mphglobal(model1, '(ez2)', 'dataset', indicator.dset,'solnum',top2Indices(1));
temp2 = mphglobal(model1, '(ez2)', 'dataset', indicator.dset,'solnum',top2Indices(2));
if temp1>temp2
    iTM = top2Indices(1);
    iTE = top2Indices(2);
else
    iTM = top2Indices(2);
    iTE = top2Indices(1);
end
check = mphglobal(model1, 'core_energy/total_energy', 'dataset', indicator.dset);
if check(iTE)>0.3
    body_TE = mphglobal(model1, '(ewfd2.beta)', 'dataset', indicator.dset,'solnum',iTE);
else
    body_TE = -1;
end
if check(iTM)>0.3
    body_TM = mphglobal(model1, '(ewfd2.beta)', 'dataset', indicator.dset,'solnum',iTM);
else
    body_TM = -1;
end

%}