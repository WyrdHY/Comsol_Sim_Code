addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Matlab Fcn Lib')));%% Initialization
%%
addpath('C:\Program Files\COMSOL\COMSOL55\Multiphysics\mli')
import com.comsol.model.*
import com.comsol.model.util.*
format long;
mphstart(2036); 
%%
modelPath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Waveguide_Vg_Matching\Straight_Wg.mph";
savePath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Waveguide_Vg_Matching\Straight_Wg_Saved.mph";

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
mesh_size = 50;
model.study('std1').feature('mode').set('neigs', '2');
disp(n_initial_guess)
%% Define some temporal function 
function loader(model,ida,n_initial_guess,core_width,core_height)
    model.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
    model.param.set('initial_guess', num2str(n_initial_guess, '%.5f'));
    model.param.set('core_width', [num2str(core_width, '%.4f'),'[um]']);
    model.param.set('core_height', [num2str(core_height, '%.4f'),'[um]']);
    if ida>800
        mesh_size = 100;
    else 
        mesh_size = 50;
    end
    model.param.set('mesh_size', [num2str(mesh_size,'%.2f'),'[nm]']);
    %disp('Configured');
end
 
function [iTE,iTM,nTE,nTM] = read_result(model,indicator)
    neff = mphglobal(model, 'real(ewfd.neff)', 'dataset', indicator.dset);
    [eigs_value, eigs] = sort(neff, 'descend');
    top2Indices = eigs(1:2);
    
    temp1 = mphglobal(model, '(ez)', 'dataset', indicator.dset,'solnum',top2Indices(1));
    temp2 = mphglobal(model, '(ez)', 'dataset', indicator.dset,'solnum',top2Indices(2));
    
    if temp1>temp2
        iTM = top2Indices(1);
        iTE = top2Indices(2);
    else
        iTM = top2Indices(2);
        iTE = top2Indices(1);
    end
    nTE = neff(iTE);
    nTM = neff(iTM);
end

function ng = group_index(x, nx)
    % Ensure x and nx are column vectors
    x = x(:);
    nx = nx(:);
    
    % Initialize the derivative array
    dndx = zeros(size(x));
    
    % Compute dn/dx using `derivest`
    for i = 1:length(x)
        dndx(i) = derivest(@(xq) interp1(x, nx, xq, 'spline'), x(i));
    end
    
    % Compute group index
    ng = nx - x .* dndx;
end
function ensureTE(result, token)
    result_key = result(token);  % Retrieve stored array
    array_diff = result_key(:,3) - result_key(:,4);  % Compute ex - ez
    % Check if any value is negative
    if any(array_diff < 0)
        warning('Potential TM mode detected at %s', token);
    end
end

loader(model,ida,n_initial_guess,core_width,core_height);
%% Run Study For TE mode
%lambda_values = [532,1064];
lambda_values = [532];
width_list = 2:0.2:3.8;
%ida_size = -10:0.5:10;
ida_size = 500:5:1100 ;
ida_size = ida_size - 532;
% Initialize dictionary (Map)
result = containers.Map();

current = 1;
total = length(width_list)*length(ida_size)*length(lambda_values);
startTime = datetime('now');  % Record start time
disp(['Start time: ', datestr(startTime, 'yyyy-mm-dd HH:MM:SS')]);

for w = width_list
    for lambda = lambda_values
        % Create the Key to hold it
        key = sprintf('%.1f_%d', w, lambda);    
        result(key) = 0;  
        ida_list = lambda+ida_size;
        
        temp = zeros(length(ida_list),4);
        temp(:,1) = ida_list; % We want n(x), this is x-array

        j = 1;
        disp([key]);
        for wavelength = ida_list
            % Load
            loader(model,wavelength,n(wavelength),w,core_height);
            model.param.set('mesh_size', [num2str(100,'%.2f'),'[nm]']);
            % Run
            indicator = flag(333); 
            studyName = indicator.std; 
            model.study(studyName).run();
            
            % Read
            [iTE,~,nTE,nTM] = read_result(model,indicator); 
            ex = mphglobal(model, '(ex)', 'dataset', indicator.dset,'solnum',iTE);
            ez = mphglobal(model, '(ez)', 'dataset', indicator.dset,'solnum',iTE);
            
            % Save
            temp(j,2) = nTE;
            temp(j,3) = ex;
            temp(j,4) = ez;
            j = j+1; 
            disp(sprintf('%d / %d',current,total));
            current = current+1;
        end
        result(key) = temp;
    end
end
endTime = datetime('now');  % Record end time
disp(['End time: ', datestr(endTime, 'yyyy-mm-dd HH:MM:SS')]);
elapsedTime = endTime - startTime;
disp(['Elapsed time: ', char(elapsedTime)]);
%%
i=1;
gap = zeros(length(width_list),1);

for w = width_list
    token_532 = sprintf('%.1f_%d',w,532);
    token_1064 = sprintf('%.1f_%d',w,1064);

    % Ensure we are comparing TE mode
    ensureTE(result, token_532);
    ensureTE(result, token_1064);
    

    % ng 532
    a_532 = result(token_532);
    x_532 = a_532(:,1);  % x values for 532 nm
    nx_532 = a_532(:,2); % Refractive index at 532 nm
    ng_532 = group_index(x_532, nx_532);
    if 0
    % ng 1064
    a_1064 = result(token_1064);
    x_1064 = a_1064(:,1);  % x values for 1064 nm
    nx_1064 = a_1064(:,2); % Refractive index at 1064 nm
    ng_1064 = group_index(x_1064, nx_1064);
    
    % Extract ng difference
    mid_i = (length(ng_532)+1)/2;
    ng_532_ = ng_532(mid_i);
    ng_1064_ = ng_1064(mid_i);
    gap(i) = ng_1064_-ng_532_;

    end
    % Plot ng_532 and ng_1064 on the same figure
    figure;
    plot(x_532,ng_1064, 'b-', 'LineWidth', 1.5); % Plot ng_1064 in blue
    hold on;
    plot(x_532,ng_532, 'r--', 'LineWidth', 1.5); % Plot ng_532 in red (dashed)
    
    % Labels and legend
    xlabel('Index in Array');
    ylabel('Group Index n_g');
    title(sprintf('Group Index for Width %.1f μm', w));
    legend('n_g at 1064 nm', 'n_g at 532 nm');
    grid on;
    hold off
    i=i+1;
end
figure;
plot(width_list,gap);
%%
keys(result)
%%
for w = width_list
    token_532 = sprintf('%.1f_%d', w, 532);
    temp = result(token_532); 
    % Extracting data
    x = temp(:,1);  % Wavelength in nm
    nx = temp(:,2); % Refractive index effective
    ex = temp(:,3); % Electric field component x
    ez = temp(:,4); % Electric field component z
    ng = group_index(x, nx); % Group Index
    figure;
    subplot(1,1,1);
    plot(x, ng, 'LineWidth', 1.5);
    xlabel('Wavelength (nm)');
    ylabel('Group Index');
    title('Group Index vs Wavelength');
    grid on;
    title(token_532);
    if 0
    subplot(3,1,2);
    plot(x, nx, 'LineWidth', 1.5);
    xlabel('Wavelength (nm)');
    ylabel('Refractive Index');
    title('Refractive Index vs Wavelength');
    grid on;
    subplot(3,1,3);
    plot(x, ex, 'r', 'LineWidth', 1.5);
    hold on;
    plot(x, ez, 'b', 'LineWidth', 1.5);
    xlabel('Wavelength (nm)');
    ylabel('Electric Field');
    title(token_532);
    legend('Ex', 'Ez');
    grid on;
    hold off;
    end
end


%%
figure; % Create a single figure for all plots
hold on; % Hold on to plot multiple curves

for w = width_list
    token_532 = sprintf('%.1f_%d', w, 532);
    temp = result(token_532);  % Ensure correct structure access

    % Extracting data
    x = temp(:,1);  % Wavelength in nm
    nx = temp(:,2); % Refractive index effective
    ng = group_index(x, nx); % Group Index

    % Plot ng vs. x with circles and dashed lines
    plot(x, ng, '-', 'LineWidth', 2, 'MarkerSize', 3, ...
         'DisplayName', sprintf('w = %.1f µm', w));
end

% Add vertical dashed lines at x = 532 and x = 1064, but remove them from legend
xline(532, '--k', '532 nm', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right', 'LineWidth', 3);
xline(1064, '--k', '1064 nm', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right', 'LineWidth', 3);

% Labeling and formatting
xlabel('Wavelength (nm)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Group Index n_g', 'FontSize', 16, 'FontWeight', 'bold');
title('Group Index vs Wavelength for Different Widths', 'FontSize', 18, 'FontWeight', 'bold');

legend('show', 'FontSize', 14, 'Location', 'best'); % Show legend with w values
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 1.5); % Increase axis font size and thickness
box on; % Ensure clear plot borders
hold off;






%%
keys(result)
%%
a = result('4.0_532');
x = a(:,1);
nx = a(:,2);
ex = a(:,3);
ez = a(:,4);
ng = group_index(x,nx);

%% Create Subplots
figure;
subplot(3,1,1); 
plot(x, ng, 'b', 'LineWidth', 1.5);
xlabel('x'); ylabel('n_x');
title('n_x vs x');
grid on;

subplot(3,1,2); 
plot(1:length(ex), ex, 'r', 'LineWidth', 1.5);
xlabel('Index in Array'); ylabel('E_x');
title('E_x vs. Index');
grid on;

subplot(3,1,3); 
plot(1:length(ez), ez, 'g', 'LineWidth', 1.5);
xlabel('Index in Array'); ylabel('E_z');
title('E_z');
grid on;

%%   
















%%
neff = mphglobal(model, 'real(ewfd.neff)', 'dataset', indicator.dset);
[~, eigs] = sort(neff, 'descend');
top2Indices = eigs(1:2);

temp1 = mphglobal(model, '(ez)', 'dataset', indicator.dset,'solnum',top2Indices(1));
temp2 = mphglobal(model, '(ez)', 'dataset', indicator.dset,'solnum',top2Indices(2));

if temp1>temp2
    iTM = top2Indices(1);
    iTE = top2Indices(2);
else
    iTM = top2Indices(2);
    iTE = top2Indices(1);
end
disp('Body Momentum Done')

modes = {'TE', 'TM'};
solnums = [iTE, iTM];  
titles = {'TE Mode, E-Right', 'TM Mode, E-Left'};
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
w0_list = (0.5:0.5:5)*1e-6;
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

