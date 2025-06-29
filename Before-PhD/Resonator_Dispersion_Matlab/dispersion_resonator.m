addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Matlab Fcn Lib')));\n%% Initialization
addpath('C:\Program Files\COMSOL\COMSOL55\Multiphysics\mli')
import com.comsol.model.*
import com.comsol.model.util.*
format long;
mphstart(2036); 

modelPath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\PG_Coupler_Jan_9_2025\Single Ring_Body.mph";
savePath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\PG_Coupler_Jan_9_2025\Single Ring_Body.mph";

model = mphload(modelPath);   
ModelUtil.showProgress(true);

%% Navigator
mphnavigator(model);
    % Std1 is for small Region ewfd2
    % Std2 is for large Region ewfd
%% Configure Parameters for Dispersion Single Ring
nClad = 1.450721;
nCore = 1.481681554; %1.498419381 = 532nm
pml_extend = 4;
width = 6;
ida = 1064;

model.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
model.param.set('nClad', num2str(nClad, '%.3f'));
model.param.set('nCore', num2str(nCore, '%.3f'));
model.param.set('pml_extend', [num2str(pml_extend, '%.2f'),'[um]']);
model.param.set('width', [num2str(width, '%.4f'),'[um]']);
%% Run Study For Dispersion
indicator = flag(0); %Do Small Region ewfd2
studyName = indicator.std; 
i=1;
model=model1;
% Scan Range
x=4;
start_wavelength = 1064-x; 
end_wavelength = 1064+x;  
step_size = 0.5;        
wavelength_range = start_wavelength:step_size:end_wavelength;

% Preallocate Array
num_points = length(wavelength_range);
wavelength = zeros(1, num_points); % Wavelengths in meters
omega_mu = zeros(1, num_points);  % Omega values
mu = zeros(1, num_points);        % Beta values
neff = zeros(1, num_points); 

for ida = wavelength_range
    model.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
    disp(['Running study: ', studyName, ' for wavelength = ', num2str(ida), ' nm']);
    model.study(studyName).run();
    disp(['             ','Done'])
    %Only Record the solnum for MaxNeff
    eigs = mphglobal(model, 'real(ewfd2.neff)', 'dataset', indicator.dset);
    [~,imax] = max(eigs);
    
    %Update the Result
    solnum = imax;
    neff(i) = eigs(imax);
    wavelength(i) = ida*1e-9;
    omega_mu(i) = mphglobal(model, 'real(ewfd2.omega)', 'dataset', indicator.dset,'solnum',solnum);
    mu(i) = mphglobal(model, 'real(ewfd2.beta)', 'dataset', indicator.dset,'solnum',solnum);
    i = i+1;
end
disp(['Study ', studyName, ' completed successfully.']);

%% Plot all of the E-Field(Should be 5)
solnum = 1;
figure;
model.result('pg1').set('data', indicator.dset);
model.result('pg1').set('solnum', num2str(solnum));
model.result('pg1').feature('surf1').set('expr', indicator.normE); % Ensure 'surf1' exists in pg1
model.result('pg1').feature('surf1').set('descr', '|E|');  
mphplot(model, 'pg1');

%%
dispersion(omega_mu,mu,1064)

%% Saving
mphsave(model, savePath);
disp(['Model saved to: ', savePath]);
%% Exit 
% ModelUtil.disconnect();
ModelUtil.remove('model');
disp('Model removed from memory.');

%%