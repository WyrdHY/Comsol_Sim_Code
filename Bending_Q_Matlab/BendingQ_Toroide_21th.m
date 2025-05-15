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
modelPath1 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Toroidal_21th\Sphere_THG.mph";
savePath1 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Toroidal_21th\Sphere_THG_save_Q.mph";


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
r_=1.5;
air_offset = 10;
pml_extend = 0.4;


%n=load_and_interpolate_n("C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\Sellmeier Fitting\Sellmeier_4%_PX332_FHD.txt");
n=load_and_interpolate_n("C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\measured_refractive_index\Ge-SiO2 (2% doping RTA) (Shijia-UCSB 2024a n 0.459-1.688 µm).txt");
nCore = n(ida);
%n = load_and_interpolate_n("C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\measured_refractive_index\ICP-PECVD Long Anneal.txt");
%n = load_and_interpolate_n("C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\measured_refractive_index\ReMeasured_Good_THOX.txt");

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
%% Calculate Bending Q. Input: ida, radius. Width is fixed.
ida_list = [1064]; 
result_map = containers.Map(); 
j = 1;
for ida = ida_list
key = sprintf("ida=%d",ida);
startTime = datetime('now');  % Record start time
disp(['Start time: ', datestr(startTime, 'yyyy-mm-dd HH:MM:SS')]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra Wavelength Layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radius_list = [2,4,6,8];

result = zeros(5,length(radius_list)); %r,TE,Q TE, TM,Q TM
result(1,:) = radius_list';
i=1;
tot = length(radius_list)*length(ida_list);
for r = radius_list
    disp(sprintf('%d/%d',j,tot));
    j=j+1;
    indicator = flag(0); % Small Region ewfd2
    studyName = indicator.std; 
    model1.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
    model1.param.set('radius', [num2str(r, '%.3f'),'[um]']);
    model1.study('std1').feature('mode').set('neigs', '6');
    model1.param.set('initial_guess', ['nCore*(radius+r_)']);
    model1.study(studyName).run();
    neff = mphglobal(model1, 'real(ewfd2.neff)', 'dataset', indicator.dset);
    core_total = mphglobal(model1, 'eigen_core_energy/eigen_total_energy', 'dataset', indicator.dset);
    % Filter out modes with core_total < 0.8
    validIndices = find(core_total >= 0.2);
    % Sort the valid eigenvalues in descending order
    [~, sortOrder] = sort(neff(validIndices), 'descend');
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
    initial_guess = neff(iTE);
    
    
    % Prepare Macro for detailed Scan
    model1.param.set('initial_guess', num2str(initial_guess));
    indicator = flag(1); % Large Region ewfd1
    studyName = indicator.std; 
    model1.study('std2').feature('mode').set('neigs', '6');
    model1.study(studyName).run();

    neff = mphglobal(model1, 'real(ewfd.neff)', 'dataset', indicator.dset);
    core_total = mphglobal(model1, 'core_energy/total_energy', 'dataset', indicator.dset);
    % Filter out modes with core_total < 0.7
    validIndices = find(core_total >= 0.2);
    watchdog = 1;
    info = sprintf("ida=%d, r=%.1f um",ida,r);

    if numel(validIndices) < 2
        disp(['Skipping ', info, ' Careful with this']);
        iTE = validIndices(1);
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

    try
        beta_TE = mphglobal(model1, 'real(ewfd.beta)', 'dataset', indicator.dset,'solnum',iTE);
        re = mphglobal(model1,'real(ewfd.neff)', 'dataset', indicator.dset,'solnum',iTE);
        im = mphglobal(model1,'imag(ewfd.neff)', 'dataset', indicator.dset,'solnum',iTE);
        tot_e_TE = mphglobal(model1, 'total_energy', 'dataset', indicator.dset, 'solnum', iTE);
        leak_e_TE = mphglobal(model1, 'leak_energy', 'dataset', indicator.dset, 'solnum', iTE);
        leak_A_TE = mphglobal(model1, 'leak_A', 'dataset', indicator.dset, 'solnum', ...
            iTE,'unit',   'um^2' );
        Q_TE = re/(im*2);
    catch
        beta_TE = 0;
        Q_TE = 0;
    end

    try 
        beta_TM = mphglobal(model1, 'real(ewfd.beta)', 'dataset', indicator.dset,'solnum',iTM);
        re = mphglobal(model1,'real(ewfd.neff)', 'dataset', indicator.dset,'solnum',iTM);
        im = mphglobal(model1,'imag(ewfd.neff)', 'dataset', indicator.dset,'solnum',iTM);
        tot_e_TM = mphglobal(model1, 'total_energy', 'dataset', indicator.dset, 'solnum', iTM);
        leak_e_TM = mphglobal(model1, 'leak_energy', 'dataset', indicator.dset, 'solnum', iTM);
        leak_A_TM = mphglobal(model1, 'leak_A', 'dataset', indicator.dset, 'solnum', ...
            iTM,'unit',   'um^2' );
        Q_TM = re/(im*2);
    catch
        beta_TM = 0;
        Q_TM = 0;
    end

    hFig=figure;
    model1.result('pg1').set('data', indicator.dset);
    model1.result('pg1').set('solnum', iTE);
    model1.result('pg1').set('edgecolor', 'white');
    model1.result('pg1').feature('surf1').set('colortable', 'Rainbow');

    log_plot = 0;
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

    model1.result('pg1').feature('arws1').active(0);
    model1.result('pg1').feature('arws1').set('expr', {'abs(ewfd.Er)','abs(ewfd.Ez)'});
    TE_field = mphplot(model1, 'pg1','rangenum',1);
    tinfo = sprintf("ida=%d, r=%.1f um,log plot=%d",ida,r,log_plot);

    title(tinfo);
    drawnow;

    % Save Temporal Photo In Case You Need It Later
    outFolder = 'C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Toroidal_21th\Temp_Output';
    filename = fullfile(outFolder, tinfo+ '.fig');
    savefig(hFig, filename);

    result(2,i) = beta_TE;  
    result(3,i) = abs(Q_TE);
    result(4,i) = leak_A_TE;  
    result(5,i) = leak_e_TE/tot_e_TE;


    result(6,i) = beta_TM;
    result(7,i) = abs(Q_TM);
    result(8,i) = leak_A_TM;
    result(9,i) = leak_e_TM/tot_e_TM;

    disp(i/length(radius_list));
    i = i+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra Wavelength Layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_map(key)=result;
end
endTime = datetime('now');  % Record end time
disp(['End time: ', datestr(endTime, 'yyyy-mm-dd HH:MM:SS')]);
elapsedTime = endTime - startTime;
disp(['Elapsed time: ', char(elapsedTime)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Bending Q for low confinement. 
% Input: ida, radius, and width
ida_list = [1064]; 
r_curve_list = [1.5,2,2.5,3];
result_map = containers.Map(); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate all of the key
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total number of combinations
nCombos = numel(ida_list)*numel(r_curve_list);

% preallocate a 1×n cell
key_holder = cell(1, nCombos);

k = 1;
for ii = 1:numel(ida_list)
    for jj = 1:numel(r_curve_list)
        ida = ida_list(ii);
        r   = r_curve_list(jj);
        key = sprintf('ida=%d_r_=%.1f', ida, r);
        % each entries{k} is a 1×3 cell: {key, ida, r}
        key_holder{k} = { key, ida, r };
        k = k + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin the loop{for each key, we loop the radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radius_list = [3.5,3.6,3.7,3.8,3.9,4,4.1,4.2,4.3];

j = 1;
startTime = datetime('now');  % Record start time
disp(['Start time: ', datestr(startTime, 'yyyy-mm-dd HH:MM:SS')]);

for key_entry = key_holder
key_entry = key_entry{1};
key = key_entry{1};
ida = key_entry{2};
r_ = key_entry{3};
    if r_ == 1.5
        radius_list = 3.5 : 0.1 : 6.5;
    elseif r_ == 2
        radius_list = 3.5 : 0.1 : 6.0;
    elseif r_ == 2.5
        radius_list = 3.5 : 0.1 : 5.5;
    elseif r_ == 3
        radius_list = 3.5 : 0.1 : 5.0;
    else
        error('Unexpected r_ = %g', r_)
    end

result = zeros(5,length(radius_list)); %r,TE,Q TE, TM,Q TM
result(1,:) = radius_list';
i=1;
tot = length(radius_list)*length(key_holder);
for r = radius_list
    disp(sprintf('%d/%d',j,tot));
    j=j+1;
    indicator = flag(0); % Small Region ewfd2
    studyName = indicator.std; 
    model1.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
    model1.param.set('radius', [num2str(r, '%.3f'),'[um]']);
    model1.param.set('r_', [num2str(r_, '%.3f'),'[um]']);
 
    model1.study('std1').feature('mode').set('neigs', '6');
    model1.param.set('initial_guess', ['nCore*(radius+r_)']);
    model1.study(studyName).run();
    neff = mphglobal(model1, 'real(ewfd2.neff)', 'dataset', indicator.dset);
    core_total = mphglobal(model1, 'eigen_core_energy/eigen_total_energy', 'dataset', indicator.dset);
    % Filter out modes with core_total < 0.8
    validIndices = find(core_total >= 0.2);
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
        continue %skip this if you cannot find a mode in easy search
    end
    initial_guess = neff(iTE);
    
    
    % Prepare Macro for detailed Scan
    model1.param.set('initial_guess', num2str(initial_guess));
    indicator = flag(1); % Large Region ewfd1
    studyName = indicator.std; 
    model1.study('std2').feature('mode').set('neigs', '6');
    model1.study(studyName).run();

    neff = mphglobal(model1, 'real(ewfd.neff)', 'dataset', indicator.dset);
    core_total = mphglobal(model1, 'core_energy/total_energy', 'dataset', indicator.dset);
    % Filter out modes with core_total < 0.7
    validIndices = find(core_total >= 0.2);
    watchdog = 1;
    info = sprintf("ida=%d, r=%.1f um",ida,r);

    if numel(validIndices) < 2
        disp(['Skipping ', info, ' Careful with this']);
        iTE = validIndices(1);
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

    try
        beta_TE = mphglobal(model1, 'real(ewfd.beta)', 'dataset', indicator.dset,'solnum',iTE);
        re = mphglobal(model1,'real(ewfd.neff)', 'dataset', indicator.dset,'solnum',iTE);
        im = mphglobal(model1,'imag(ewfd.neff)', 'dataset', indicator.dset,'solnum',iTE);
        tot_e_TE = mphglobal(model1, 'total_energy', 'dataset', indicator.dset, 'solnum', iTE);
        leak_e_TE = mphglobal(model1, 'leak_energy', 'dataset', indicator.dset, 'solnum', iTE);
        leak_A_TE = mphglobal(model1, 'leak_A', 'dataset', indicator.dset, 'solnum', ...
            iTE,'unit',   'um^2' );
        Q_TE = re/(im*2);
    catch
        beta_TE = 0;
        Q_TE = 0;
    end

    try 
        beta_TM = mphglobal(model1, 'real(ewfd.beta)', 'dataset', indicator.dset,'solnum',iTM);
        re = mphglobal(model1,'real(ewfd.neff)', 'dataset', indicator.dset,'solnum',iTM);
        im = mphglobal(model1,'imag(ewfd.neff)', 'dataset', indicator.dset,'solnum',iTM);
        tot_e_TM = mphglobal(model1, 'total_energy', 'dataset', indicator.dset, 'solnum', iTM);
        leak_e_TM = mphglobal(model1, 'leak_energy', 'dataset', indicator.dset, 'solnum', iTM);
        leak_A_TM = mphglobal(model1, 'leak_A', 'dataset', indicator.dset, 'solnum', ...
            iTM,'unit',   'um^2' );
        Q_TM = re/(im*2);
    catch
        beta_TM = 0;
        Q_TM = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hFig=figure;
    model1.result('pg1').set('data', indicator.dset);
    model1.result('pg1').set('solnum', iTE);
    model1.result('pg1').set('edgecolor', 'white');
    model1.result('pg1').feature('surf1').set('colortable', 'Rainbow');

    log_plot = 1;
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

    model1.result('pg1').feature('arws1').active(0);
    model1.result('pg1').feature('arws1').set('expr', {'abs(ewfd.Er)','abs(ewfd.Ez)'});
    TE_field = mphplot(model1, 'pg1','rangenum',1);
    tinfo = sprintf("%s,Radius=%.1f um,log plot=%d",key,r,log_plot);

    title(tinfo, 'Interpreter', 'none')
    drawnow;

    % Save Temporal Photo In Case You Need It Later
    outFolder = 'C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Toroidal_21th\Temp_Output';
    filename = fullfile(outFolder, tinfo+ '.fig');
    savefig(hFig, filename);
    close;
    result(2,i) = beta_TE;  
    result(3,i) = abs(Q_TE);
    result(4,i) = leak_A_TE;  
    result(5,i) = leak_e_TE/tot_e_TE;


    result(6,i) = beta_TM;
    result(7,i) = abs(Q_TM);
    result(8,i) = leak_A_TM;
    result(9,i) = leak_e_TM/tot_e_TM;

    disp(i/length(radius_list));
    i = i+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra Wavelength Layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_map(key)=result;
end
endTime = datetime('now');  % Record end time
disp(['End time: ', datestr(endTime, 'yyyy-mm-dd HH:MM:SS')]);
elapsedTime = endTime - startTime;
disp(['Elapsed time: ', char(elapsedTime)]);
%%
protect = result_map;
%%
a = protect('ida=1064');
%% Plot Bending Q at specific key
ida = 1064;
r_ = 1.5;
key = sprintf('ida=%d_r_=%.1f', ida, r_);
%key = sprintf('ida=%d', ida);
result = a %result_map(key);
r      = result(1,:);   % r values
betaTE = result(2,:);
Q_TE   = abs(result(3,:));
betaTM = result(4,:);
Q_TM   = abs(result(5,:));
figure;
loglog(r, Q_TE, '--d', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
hold on;
yline(1e9, '--', 'LineWidth', 1);grid on;
hold off;
xlabel('r');
ylabel('Q');
title('21th Toroid, 450nm');

%% Plot all of them 
figure;  
for key = key_holder
    key = key{1};
    key = key{1};
    try
    result  = result_map(key);             
    r       = result(1,:);                
    Q_TE    = abs(result(3,:));           
    loglog(r, Q_TE, '-o', ...
           'LineWidth', 1.5, ...
           'DisplayName', key);hold on;
    end
end

yline(1e9, '--', 'LineWidth', 1, 'DisplayName', 'Q = 10^9');
%ylim([1e4,3e9]);
xlabel('r (µm)');
ylabel('Q');
title('r vs. Q (log–log scale)');
legend('show', 'Interpreter', 'none');
grid on;  hold off;
%% Save the Result
out = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Toroidal_21th\Temp_Output";

%--- Step 1: Pull your results ---------------------------------------------
key = 'ida=450';
result1 = protect(key);


%--- Step 2: Convert to tables --------------------------------------------
makeTable = @(res) table( ...
    res(1,:)', ...                       % radius
    res(2,:)', ...                       % beta_TE
    res(3,:)', ...                       % Q_TE
    res(4,:)', ...                       % leak_A_TE
    res(5,:)', ...                       % leak_e_TE ./ tot_e_TE
    'VariableNames', { ...
        'radius', ...
        'beta_TE', ...
        'Q_TE(unit=1)', ...
        'leak_A_TE(um^2)', ...
        'leak_ratio_TE' ...
    } ...
);

T1 = makeTable(result1);


%--- Step 3: Write to one Excel workbook, separate sheets -----------------
filename = 'all_ida_results.xlsx';
output = fullfile(out, filename);
writetable(T1, output, 'Sheet', key);

fprintf('Wrote tables to %s \n', output);




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