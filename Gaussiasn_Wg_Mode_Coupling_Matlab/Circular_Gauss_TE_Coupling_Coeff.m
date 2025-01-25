addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Matlab Fcn Lib')));%% Initialization
which Fcn_Lib.m;
%%
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Matlab Fcn Lib')));
disp(path);  % Display current MATLAB path to verify
%%
if exist('load_and_interpolate_n', 'file')
    disp('Function load_and_interpolate_n is accessible.');
else
    disp('Function load_and_interpolate_n is NOT accessible.');
end

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
mphnavigator(model);

%% Configure Parameters for Dispersion Single Ring
fpath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Material\Sellmeier Fitting\Sellmeier_2%_Long Anneal.txt";
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

