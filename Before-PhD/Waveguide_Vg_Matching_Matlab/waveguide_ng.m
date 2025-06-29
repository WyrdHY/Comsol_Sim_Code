%% Include the Function Lib
basefolder = "D:\Caltech\Comsol_Simulation\Code\Matlab Fcn Lib";
addpath(genpath(basefolder)); 
addpath('D:\COMSOL 6_2\COMSOL62\Multiphysics\mli')
import com.comsol.model.*
import com.comsol.model.util.*
format long;
mphstart(2036); 
%%
modelPath = "D:\Caltech\Comsol_Simulation\Model\Straight_Wg.mph";
savePath = "D:\Caltech\Comsol_Simulation\Model\Straight_Wg_ng.mph";

model = mphload(modelPath);   
com.comsol.model.util.ModelUtil.showProgress(true);

%% Navigator
mphnavigator(model);
%% Run Study For TE mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda_ = 1550;
lambda_mask = -1:1:1;
lambda_values = lambda_+lambda_mask;
lambda_values = [980,1550];
model.param.set('width', [num2str(10, '%.3f'),'[um]']);
model.param.set('height', [num2str(1.6, '%.3f'),'[um]']);
model.param.set('mesh_size', [num2str(100, '%.3f'),'[nm]']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The result is a tables 3xn
% ida, nTE, nTM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = zeros(length(lambda_values),3);
result(:,1) = lambda_values;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure Loop Para
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current = 1;
total = length(lambda_values);
startTime = datetime('now');  % Record start time
disp(['Start time: ', datestr(startTime, 'yyyy-mm-dd HH:MM:SS')]);

indi = struct();
indi.dset = 'dset1';
indi.std = 'std1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
j = 1;
for lambda = lambda_values

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model.param.set('ida', [num2str(lambda, '%.3f'),'[nm]']);
    key = sprintf("ida=%.3fnm",lambda);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [iTE,iTM,skip] = calculate(model);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read and Save
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nTE = mphglobal(model, 'real(ewfd.neff)', 'dataset', 'dset1','solnum',iTE);
    nTM = mphglobal(model, 'real(ewfd.neff)', 'dataset', 'dset1','solnum',iTM);

    result(j,2) = nTE;
    result(j,3) = nTM;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    indi = struct();
    indi.dset = 'dset1';
    indi.std = 'std1';
    log_plot = 0;
    tinfo = ['TE ',key];
    arrow =0;
    figure;
    visualize(model,indi,iTE,log_plot,arrow,tinfo);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Next Run
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    j = j+1; 
    disp(sprintf('%d / %d',current,total));
    current = current+1;
end

endTime = datetime('now');  % Record end time
disp(['End time: ', datestr(endTime, 'yyyy-mm-dd HH:MM:SS')]);
elapsedTime = endTime - startTime;
disp(['Elapsed time: ', char(elapsedTime)]);

%%
a = result;
x = a(:,1);          
nTE = a(:,2);       
nTM = a(:,3);        


ng_TE = group_index(x, nTE);
ng_TM = group_index(x, nTM);

figure;
subplot(2,1,1);
plot(x, nTE, 'b-', 'LineWidth', 1.5); hold on;
plot(x, nTM, 'r--', 'LineWidth', 1.5);
xlabel('Wavelength (nm)', 'FontSize', 12);
ylabel('Refractive Index (n)', 'FontSize', 12);
title('Effective Refractive Index vs Wavelength', 'FontSize', 14);
legend('n_{TE}', 'n_{TM}', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(x, ng_TE, 'b-', 'LineWidth', 1.5); hold on;
plot(x, ng_TM, 'r--', 'LineWidth', 1.5);
xlabel('Wavelength (nm)', 'FontSize', 12);
ylabel('Group Index (n_g)', 'FontSize', 12);
title('Group Index vs Wavelength', 'FontSize', 14);
legend('n_{g,TE}', 'n_{g,TM}', 'Location', 'best');
grid on;

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











%% Saving
mphsave(model, savePath);
disp(['Model saved to: ', savePath]);
%% Exit 
% ModelUtil.disconnect();
ModelUtil.remove('model');
disp('Model removed from memory.');

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

function [iTE,iTM,skip]=calculate(model1) 
    % As for wg, there is no need for two trials
    % I didn't use indicator but directly call the tag
    skip = 0;
    model1.param.set('initial_guess', 'nCore*1.00');
    model1.study('std1').feature('mode').set('neigs', '2');
    model1.study('std1').run();

    neff = mphglobal(model1, 'real(ewfd.neff)', 'dataset', 'dset1');
    core_total = mphglobal(model1, 'core_energy/total_energy', 'dataset', 'dset1');
    validIndices = find(core_total >= 0.2);
    watchdog = 1;

    if numel(validIndices) < 1
        disp(['Skipping ' ' Careful with this']);
        iTE = validIndices(1);
        iTM = 0;
        xxx = mphglobal(model1, '(ex)', 'dataset', 'dset1','solnum',iTE);
        watchdog = 0;
        skip = 1;
    end
    
    if watchdog
        % Sort the valid eigenvalues in descending order
        [~, sortOrder] = sort(neff(validIndices), 'descend');
        top2Indices = validIndices(sortOrder(1:2));
    
        temp1 = mphglobal(model1, '(ez)', 'dataset', 'dset1','solnum',top2Indices(1));
        temp2 = mphglobal(model1, '(ez)', 'dataset', 'dset1','solnum',top2Indices(2));
        if temp1>temp2
            iTM = top2Indices(1);
            iTE = top2Indices(2);
        else
            iTM = top2Indices(2);
            iTE = top2Indices(1);
        end
    end

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
    model1.result('pg1').feature('arws1').set('expr', {'abs(ewfd.Ex)','abs(ewfd.Ey)'});
    TE_field = mphplot(model1, 'pg1','rangenum',1);

    title(tinfo, 'Interpreter', 'none')
    drawnow;
end