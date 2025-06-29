%% Initialize
% Model Path
opt1 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\BendQ_2025_06_25\One Ring.mph";
opt2 = "D:\Caltech\Comsol_Simulation\Model\One Ring.mph";
% Fuc Path
b1 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Matlab Fcn Lib\Fcn_lib";
b2 = "D:\Caltech\Comsol_Simulation\Code\Matlab Fcn Lib";
% Comsol Path
c1 = 'C:\Program Files\COMSOL\COMSOL62\Multiphysics\mli';
c2 = 'D:\COMSOL 6_2\COMSOL62\Multiphysics\mli';
[model_path, save_path, basepath, comsol_path] = choosePaths(opt1, opt2, b1, b2, c1, c2);

basefolder = basepath;
addpath(genpath(basefolder)); 

% Clean up workspace
clear opt1 opt2 b1 b2 c1 c2
disp(model_path)
disp(save_path)
disp(basepath)
disp(comsol_path)

%% Initialization
addpath(comsol_path)
import com.comsol.model.*
import com.comsol.model.util.*
format long;
try
mphstart(2036); 
end

modelPath1 =model_path;
savePath1 = save_path;
model1 = ModelUtil.create('Model1');
model1 = mphload(modelPath1, 'Model1');
ModelUtil.showProgress(true);
disp('Done')
%% Navigator
mphnavigator(model1);
% Std1 is for small Region ewfd2
    % Std2 is for large Region ewfd
%%


%% Global Parameters Test
ida=1550;
radius=1000;
air_offset = 30;
r=200;
width = 10;
pml_extend = 3;
height = 4; 
meshsize = 100;
air_left = 5;
thox  = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\measured_refractive_index\ReMeasured_Good_THOX.txt";
rta2 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\measured_refractive_index\Ge-SiO2 (2% doping RTA) (Shijia-UCSB 2024a n 0.459-1.688 µm).txt";
long2 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\measured_refractive_index\ICP-PECVD Long Anneal.txt";
fhd4 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\measured_refractive_index\PX332-3(4% FHD).txt";

n=load_and_interpolate_n(thox);
nClad = n(ida);
n = load_and_interpolate_n(long2);
nCore = n(ida);
betaR(nClad, nCore, ida, r, air_offset, width, 1)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ida=1550;
radius=1500;
air_offset = 20;
width = 12;
height = 2; 

pml_extend = 4;
s_scale = 1;
s_curve = 1;
meshsize = 100;

model1.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
model1.param.set('radius', [num2str(radius, '%.3f'),'[um]']);
model1.param.set('core_width', [num2str(width, '%.4f'),'[um]']);
model1.param.set('core_height', [num2str(height, '%.4f'),'[um]']);
model1.param.set('air_offset', [num2str(air_offset, '%.4f'),'[um]']);
model1.param.set('pml_extend', [num2str(pml_extend, '%.2f'),'[um]']);

model1.param.set('core_mesh_size', [num2str(meshsize),'[nm]']);
model1.study('std2').feature('mode').set('eigwhich', 'lm');  
model1.param.set('s_scale', num2str(s_scale)); % set s
model1.param.set('s_curve', num2str(s_curve)); % set s_curve
disp('Para Configured')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametric Sweep List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale_and_curve = [
    [1,1];
    ];
THOX_list = [5];
Si_list = [1];
radius_list = [500];
j = 1;
tot = length(radius_list)*length(Si_list)*length(THOX_list)*size(scale_and_curve,1);
result_total = containers.Map();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startTime = datetime('now');  % Record start time
disp(['Start time: ', datestr(startTime, 'yyyy-mm-dd HH:MM:SS')]);

for idx_pml = 1:size(scale_and_curve,1)
    s_val = scale_and_curve(idx_pml, 1);
    s_curve_val = scale_and_curve(idx_pml, 2);
    
    model1.param.set('s_scale', num2str(s_val)); % set s
    model1.param.set('s_curve', num2str(s_curve_val)); % set s_curve
    for Si = Si_list
        model1.material('mat17').active(Si);
        for THOX = THOX_list
        key = sprintf("ida=%.1f, THOX=%.1fum, Si=%d, s=%d, s_curve=%d", ida, THOX, Si, s_val, s_curve_val); disp(key);
        model1.param.set('THOX_height', [num2str(THOX, '%.4f'),'[um]']);
        result = zeros(3,length(radius_list)); %r,loss,Q_TE, 
        result(1,:) = radius_list';
        i = 1;
            for r = radius_list
                %%%%%%%%%%%%%%%%%%
                % Run 
                %%%%%%%%%%%%%%%%%%
                model1.param.set('radius', [num2str(r, '%.3f'),'[um]']);
                iTE = calculate(model1);
                %%%%%%%%%%%%%%%%%%
                % Read
                %%%%%%%%%%%%%%%%%%
                beta_TE = mphglobal(model1, 'real(ewfd.beta)', 'dataset', 'dset2','solnum',iTE);
                re = mphglobal(model1,'real(ewfd.neff)', 'dataset', 'dset2','solnum',iTE);
                im = mphglobal(model1,'imag(ewfd.neff)', 'dataset', 'dset2','solnum',iTE);
                loss_db_m = mphglobal(model1,'20/log(10)*imag(ewfd.neff)*ewfd.k0', 'dataset', 'dset2','solnum',iTE);
                Q_TE = abs(re/(im*2));  
                
                if 1
                    % --- Read power monitors ---
                    BndPrPML = mphglobal(model1, 'BndPrPML', 'dataset', 'dset2', 'solnum', iTE);
                    DomPrPML = mphglobal(model1, 'DomPrPML', 'dataset', 'dset2', 'solnum', iTE);
                    Bnd0     = mphglobal(model1, 'Bnd0',     'dataset', 'dset2', 'solnum', iTE);
                    Dom0     = mphglobal(model1, 'Dom0',     'dataset', 'dset2', 'solnum', iTE);
                
                    % --- Calculate absorption efficiencies (%) ---
                    eps_val   = 1e-30;  % avoid division by zero
                    effBnd    = (Bnd0   - BndPrPML) ./ max(Bnd0, eps_val) * 100;
                    effDom    = (Dom0   - DomPrPML) ./ max(Dom0, eps_val) * 100;
                
                    % --- Print all power monitor values with labels + efficiencies ---
                    fprintf('\nPower Monitoring Results (solnum = %d)\n', iTE);
                    fprintf('--------------------------------------------------\n');
                    fprintf('  ▸ Bnd0       (before PML, boundary)   : %.4e  W\n', Bnd0);
                    fprintf('  ▸ BndPrPML   (inside  PML, boundary) : %.4e  W\n', BndPrPML);
                    fprintf('    → Absorption Efficiency (boundary): %.2f %%\n\n', effBnd);
                
                    fprintf('  ▸ Dom0       (before PML, domain)     : %.4e  W\n', Dom0);
                    fprintf('  ▸ DomPrPML   (inside  PML, domain)   : %.4e  W\n', DomPrPML);
                    fprintf('    → Absorption Efficiency (domain)  : %.2f %%\n', effDom);
                    fprintf('--------------------------------------------------\n\n');
                end

                %%%%%%%%%%%%%%%%%%
                % Plot 
                %%%%%%%%%%%%%%%%%%
                log_plot = 1;
                arrow = 1;
                tinfo = sprintf("radius = %d, %s", r, key);
                figure;
                visualize(model1,flag(1),iTE,log_plot,arrow,tinfo)
                %%%%%%%%%%%%%%%%%%
                % Save 
                %%%%%%%%%%%%%%%%%%
                result(2,i) = loss_db_m;
                result(3,i) = Q_TE;
                %%%%%%%%%%%%%%%%%%
                % Next 
                %%%%%%%%%%%%%%%%%%
                i = i+1;
                j = j+1;fprintf('\n%d/%d\n',j,tot);
            end
        result_total(key) = result;
        end
    end
end

mphsave(model1, savePath1);
endTime = datetime('now');  % Record end time
disp(['End time: ', datestr(endTime, 'yyyy-mm-dd HH:MM:SS')]);
elapsedTime = endTime - startTime;
disp(['Elapsed time: ', char(elapsedTime)]);
disp(repmat('-', 1, 40));

%% systematic plot
fs = 14;                       

figure;
set(gca,'FontSize',fs);

% --- Extract & sort keys by THOX descending ---
keys   = result_total.keys;      
nK     = numel(keys);
THOX_v = zeros(1, nK);
for k = 1:nK
    t = regexp(keys{k}, 'THOX=(\d+\.?\d*)um', 'tokens', 'once');
    THOX_v(k) = str2double(t{1});
end
[~, sortIdx] = sort(THOX_v, 'descend');
sortedKeys   = keys(sortIdx);

% --- Parse out the constants from the first key ---
firstKey = sortedKeys{1};
idaVal      = str2double(regexp(firstKey, 'ida=([\d\.]+)',     'tokens','once'));
sVal        = str2double(regexp(firstKey, 's=(\d+)',            'tokens','once'));
sCurveVal   = str2double(regexp(firstKey, 's_curve=(\d+)',      'tokens','once'));
SiVal_const = str2double(regexp(firstKey, 'Si=(\d+)',           'tokens','once'));

base2col = containers.Map;     % map base → colour index
colours   = lines(nK);         
markers   = {'o','s','d','^','v','>','<','p','h'};

for idx = 1:nK
    key  = sortedKeys{idx};
    dat  = result_total(key);      % 3×N  [r; loss; Q_TE]

    r     = dat(1,:);
    Q_TE  = abs(dat(3,:));

    % Extract only THOX for the legend label
    thoxVal = regexp(key, 'THOX=(\d+\.?\d*)um', 'tokens','once');
    thoxVal = str2double(thoxVal{1});

    base = sprintf('THOX = %.1f µm', thoxVal);

    if ~isKey(base2col, base)
        base2col(base) = length(base2col)+1;
    end
    cIdx = base2col(base);
    col  = colours(cIdx,:);

    lstyle = '-';              % assume Si=1; if you had Si=0 variants, you could still dashed here
    mIdx   = mod(cIdx-1, numel(markers)) + 1;
    mk     = markers{mIdx};

    semilogy(r, Q_TE, lstyle, ...
        'Color',           col, ...
        'Marker',          mk, ...
        'MarkerFaceColor', 'none', ...
        'LineWidth',       1.5, ...
        'DisplayName',     base);
    hold on;
end

% Reference line
yline(1e9, '--k', 'LineWidth', 1, 'DisplayName', 'Q = 10^9');

xlabel('r (µm)', 'FontSize', fs);
ylabel('Q_{TE}', 'FontSize', fs);

% Put the constants into the title, after a newline:
title(sprintf('r vs Q_{TE} (Log Scale)\nida = %.1f nm, s = %d, s_curve = %d, Si = %d', ...
    idaVal, sVal, sCurveVal, SiVal_const), ...
    'FontSize', fs+2);

legend('Location','best','FontSize',fs);
grid on;

%% Visualize it
log_plot = 0;
arrow = 1;
tinfo = "1";
figure;
visualize(model1,flag(0),iTE,log_plot,arrow,tinfo)
%%
indicator = flag(1);
ex = mphglobal(model1, 'ex', 'dataset', indicator.dset);
    
ez = mphglobal(model1, 'ez', 'dataset', flag(1).dset);
core_total = mphglobal(model1, 'core_energy/total_energy', 'dataset', indicator.dset);

% ——— 2) Filter 1: core_total ≥ 0.2 ———
valid1  = find(core_total >= 0.001);

% ——— 3) Filter 2: ex/ez ≥ 1 ———
ratios = ex(valid1) ./ ez(valid1);
valid2 = valid1(ratios >= 1);


% ——— 4) Select top-2 by ex descending ———
[~, sortOrder] = sort(ez(valid2), 'descend');
top2 = valid2(sortOrder(1:min(2,numel(sortOrder))));
%% Saving
if 1

    mphsave(model1, savePath1);
    disp(['Model saved to: ', savePath1]);

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
function [model_path, save_path, basepath, comsol_path] = choosePaths(opt1, opt2, b1, b2, c1, c2)
    % Choose model_path
    if isfile(opt1)
        model_path = opt1;
    else
        model_path = opt2;
    end

    % Derive save_path
    [folder, name, ext] = fileparts(model_path);
    save_path = fullfile(folder, name + "_saved" + ext);

    % Choose basepath
    if isfolder(b1)
        basepath = b1;
    else
        basepath = b2;
    end

    % Choose COMSOL installation path
    if isfolder(c1)
        comsol_path = c1;
    else
        comsol_path = c2;
    end
end
function [iTE,iTM,skip]=eigen_locate(model1,indicator)
    skip = 0;
    neff = mphglobal(model1, 'real(ewfd2.neff)', 'dataset', indicator.dset);
    core_total = mphglobal(model1, 'eigen_core_energy/eigen_total_energy', 'dataset', indicator.dset);
    % Filter out modes with core_total < 0.6
    validIndices = find(core_total >= 0.08);
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
        [~,iTE]=max(neff);
        iTM=0;
        skip = 1; 
    end
end
function [iTE]=calculate(model1)
    % Prepare Macro for Simulation 
    model1.param.set('initial_guess', ['nCore*1.001']);
    indicator = flag(0); % Small Region ewfd2
    studyName = indicator.std; 
    model1.study('std1').feature('mode').set('neigs', '4');
    model1.study(studyName).run();
    [iTE,iTM,skip] = eigen_locate(model1,indicator);

    
    neff = mphglobal(model1, 'real(ewfd2.neff)', 'dataset', indicator.dset);
    initial_guess = neff(iTE);
    % Prepare Macro for detailed Scan
    model1.param.set('initial_guess', num2str(initial_guess));
    indicator = flag(1); % Large Region ewfd1
    studyName = indicator.std; 
    model1.study('std2').feature('mode').set('neigs', '2');
    model1.study(studyName).run();

    neff = mphglobal(model1, 'real(ewfd.neff)', 'dataset', indicator.dset);
    core_total = mphglobal(model1, 'core_energy/total_energy', 'dataset', indicator.dset);
    ex = mphglobal(model1, 'ex', 'dataset', indicator.dset);
    ez = mphglobal(model1, 'ez', 'dataset', indicator.dset);

    % ——— 2) Filter 1: core_total ≥ 0.2 ———
    valid1  = find(core_total >= 0.02);
    
    % ——— 3) Filter 2: ex/ez ≥ 1 ———
    ratios = ex(valid1) ./ ez(valid1);
    valid2 = valid1(ratios >= 1);


    % ——— 4) Select top-2 by ex descending ———
    [~, sortOrder] = sort(ez(valid2), 'descend');
    top2 = valid2(sortOrder(1:min(2,numel(sortOrder))));
    % Assign your mode indices:
    if numel(top2) >= 2
        iTM = top2(1);
        iTE = top2(2);
    elseif numel(top2) == 1
        iTE = top2(1);
        iTE2 = 0;    % only one valid mode
    else
        [~,iTE]=max(core_total);
        iTE2 = 0;
    end
    %fprintf("iTE=%d, neff=%.3f",iTE,initial_guess);

    %check = mphglobal(model1, 'core_energy/total_energy', 'dataset', indicator.dset);
    
    %[~,iTE] = max(check);
end