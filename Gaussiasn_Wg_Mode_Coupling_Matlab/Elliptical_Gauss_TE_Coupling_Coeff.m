addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Matlab Fcn Lib')));%% Initialization
%%
addpath('C:\Program Files\COMSOL\COMSOL55\Multiphysics\mli')
import com.comsol.model.*
import com.comsol.model.util.*
format long;
mphstart(2036); 
%%
modelPath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Lensed_Fiber_to_Waveguide_Coupling\Straight_Wg.mph";
savePath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Lensed_Fiber_to_Waveguide_Coupling\Straight_Wg.mph";
model = mphload(modelPath);   
ModelUtil.showProgress(true);

%% Navigator
%mphnavigator(model);

% This simulation is for
% PECVD5, Cladding, 532nm, 10*2 um
%% Configure Parameters for Dispersion Single Ring
fpath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Material\Sellmeier Fitting\Sellmeier_2%_Long Anneal.txt";
n = load_and_interpolate_n(fpath);
core_width = 10;
core_height = 2;
ida = 520;
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
z_list = (100:10:400);
theta_x=22;
theta_y=8;



w0_list = (0.5:0.2:5)*1e-6;
polar_list = [1];
eta_list = zeros(1,length(z_list));

for polar = polar_list
    i=1;
    for z = z_list
        % 0-iTM, 1-iTE
        if polar ==1 
            solnum=iTE;
            j = 1;
        else
            solnum=iTM;
            j=2;
        end
        %Para for Gaussian Beam
        nx = 1000;  
        ny = 800;  
        %Meshing 
        rx = 14.9*1e-6;
        ry = 8*1e-6;
        xVec = linspace(-rx, rx, nx);
        yVec = linspace(-ry, ry, ny);
        [X, Y] = meshgrid(xVec, yVec); 
        coords = [X(:)'; Y(:)'];
        %Extract
        Ex_data = reshape(mphinterp(model, 'ewfd.Ex','coord',coords,'dataset',indicator.dset,'solnum',solnum), ny, nx);
        Ey_data = reshape(mphinterp(model, 'ewfd.Ey','coord',coords,'dataset',indicator.dset,'solnum',solnum), ny, nx);
        %Generate Gaussian 
        [E_gauss_x, E_gauss_y] = Ellptical_Gauss_z(X,Y,polar,z,theta_x,theta_y,ida); 
        % Calculate spatial step sizes
        dx = xVec(2) - xVec(1);
        dy = yVec(2) - yVec(1);
        dA = dx * dy;  % Elemental area
        % Compute the numerator: |∫ E1* E2 dA|^2
        I12 = sum(conj(Ex_data).*E_gauss_x + conj(Ey_data).*E_gauss_y, 'all') * dA;
        numerator = abs(I12)^2;
        
        % Compute the denominators: ∫ |E1|^2 dA and ∫ |E2|^2 dA
        I11 = sum(abs(Ex_data).^2 + abs(Ey_data).^2, 'all') * dA;
        I22 = sum(abs(E_gauss_x).^2 + abs(E_gauss_y).^2, 'all') * dA;
        
        % Compute mode overlap efficiency η
        eta_list(1, i) = numerator / (I11 * I22);
        
        % Increment index for storing results
        i = i + 1;
    end
end
%% Plot ita v.s. z
te_color = [0 0.447 0.741];   % Nice blue color for TE
tm_color = [0.850 0.325 0.098]; % Nice red-orange color for TM

figure;
hold on;

plot(z_list, eta_list, '--d', ...
    'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', te_color, 'Color', te_color, 'DisplayName', 'TE Mode');


xlabel('z(gap distance) [\mum]', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Coupling Efficiency \eta', 'FontSize', 16, 'FontWeight', 'bold');
title('Coupling Efficiency vs z', 'FontSize', 18, 'FontWeight', 'bold');

grid on;
set(gca, 'GridLineStyle', '--', 'FontSize', 14);

legend('Location', 'northwest', 'FontSize', 14);

hold off;
%%




%% Debug

% Parameters for Gaussian Beam
nx = 800;  
ny = 300;  
% Meshing
rx = 14.9*1e-6;
ry = 8*1e-6;
xVec = linspace(-rx, rx, nx);
yVec = linspace(-ry, ry, ny);
[X, Y] = meshgrid(xVec, yVec); 
coords = [X(:)'; Y(:)'];
% Extract COMSOL fields
Ex_data = reshape(mphinterp(model, 'ewfd.Ex', 'coord', coords, 'dataset', indicator.dset, 'solnum', solnum), ny, nx);
Ey_data = reshape(mphinterp(model, 'ewfd.Ey', 'coord', coords, 'dataset', indicator.dset, 'solnum', solnum), ny, nx);

%Compute
polarization = 1; z = 200; theta_x = 22;theta_y=8;ida=532;
Gauss_polar=polarization;
[E_gauss_x, E_gauss_y] = Ellptical_Gauss_z(X,Y,polarization,z,theta_x,theta_y,ida); 

%
xVec_um = xVec*1e6;
yVec_um = yVec*1e6;

% Plot Gaussian field
% Set dynamic title based on polarization
figure;
subplot(1,2,1);
imagesc(xVec_um, yVec_um, sqrt(E_gauss_y.^2+E_gauss_x.^2));

set(gca, 'YDir', 'normal');  
axis image;                 
colorbar;
xlabel('x [\mum]', 'FontSize', 12);
ylabel('y [\mum]', 'FontSize', 12);
title(sprintf('laser diode, z=%.2f um',z))
% Plot COMSOL field
subplot(1,2,2);
if Gauss_polar == 1
    imagesc(xVec_um, yVec_um, abs(Ey_data));
else 
    imagesc(xVec_um, yVec_um, abs(Ex_data));
end
set(gca, 'YDir', 'normal');  
title('COMSOL Simulation')
xlabel('x [\mum]', 'FontSize', 12);
ylabel('y [\mum]', 'FontSize', 12);
axis image;                  
colorbar;



%% Saving
mphsave(model, savePath);
disp(['Model saved to: ', savePath]);
%% Exit 
% ModelUtil.disconnect();
ModelUtil.remove('model');
disp('Model removed from memory.');

