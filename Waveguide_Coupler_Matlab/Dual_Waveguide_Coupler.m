addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Matlab Fcn Lib')));%% Initialization
%%
addpath('C:\Program Files\COMSOL\COMSOL55\Multiphysics\mli')
import com.comsol.model.*
import com.comsol.model.util.*
format long;
mphstart(2036); 
%%
modelPath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Two_Waveguide_Coupler\waveguide_coupler.mph";
savePath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Two_Waveguide_Coupler\waveguide_coupler_saved.mph";

model = mphload(modelPath);   
ModelUtil.showProgress(true);

%% Navigator
mphnavigator(model);

%% Configure Parameters for Two Waveguides
fpath = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\Material\Sellmeier Fitting\Sellmeier_2%_Long Anneal.txt";

n = load_and_interpolate_n(fpath);
core_width = 9;
core_height = 2;
ida = 1550;
disp(n(ida));
n_initial_guess = 1.471; %n(ida);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trace how kappa varies with width
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  key = sprintf('w_%.1f_h_%.1f_%.1f', w, h,delta);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare Looping Element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
width_list = [10];
height_list = [4];
delta_list = [5,7,10,11,15,16,17,18];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize dictionary (Map)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = containers.Map();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
model.param.set('initial_guess', num2str(n_initial_guess, '%.5f'));
model.param.set('mesh_size', '50[nm]');
model.param.set('mesh_size2', '50[nm]');
model.study('std1').feature('mode').set('neigs', '4');
PECVD = 0;
if PECVD
    model.material('mat1').active(0)
    model.material('mat4').active(1)
else
    model.material('mat1').active(1)
    model.material('mat4').active(0)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Record Computation Time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current = 1;
total = length(width_list)*length(height_list)*length(delta_list);
startTime = datetime('now');  % Record start time
disp(['Start time: ', datestr(startTime, 'yyyy-mm-dd HH:MM:SS')]);

for w = width_list
    for h = height_list

        % Create the Key to hold it
        key = sprintf('w_%.1f_h_%.1f', w, h);    
        result(key) = 0;  

        % We want kappa(delta)
        temp = zeros(length(delta_list),1);

        j = 1;
        disp(key);
        for delta = delta_list
            disp(sprintf('%d / %d',current,total));
            % Load
            model.param.set('core_width', [num2str(w, '%.3f'),'[um]']);
            model.param.set('core2_width', [num2str(w, '%.3f'),'[um]']);
            model.param.set('core_height', [num2str(h, '%.3f'),'[um]']);
            model.param.set('delta', [num2str(delta, '%.3f'),'[um]']);


            % Run
            indicator = flag(333); 
            studyName = indicator.std; 
            model.study(studyName).run();
            
            % Read
            beta = mphglobal(model, '(ewfd.beta)', 'dataset', indicator.dset);
            ex = mphglobal(model, '(ex)', 'dataset', indicator.dset);
            [~, idx] = sort(ex, 'descend');
            n1 = idx(1);
            n2 = idx(2);
            kappa = abs(beta(n1) - beta(n2));


            % Save
            temp(j,1) = kappa;
            j = j+1; 
            current = current+1;
        end
        result(key) = temp;
    end
end
endTime = datetime('now');  % Record end time
disp(['End time: ', datestr(endTime, 'yyyy-mm-dd HH:MM:SS')]);
elapsedTime = endTime - startTime;
disp(['Elapsed time: ', char(elapsedTime)]);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject = result;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = height_list(1);

figure;
hold on;
legend_entries = {};

for w = width_list
    key = sprintf('w_%.1f_h_%.1f', w, h); 
    y = subject(key);
    plot(delta_list, y, '--d', 'LineWidth', 2, 'MarkerSize', 8);
    legend_entries{end+1} = sprintf('w = %.1f', w);
end

xlabel('$\delta$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\kappa$', 'Interpreter', 'latex', 'FontSize', 16);
yscale log;
legend(legend_entries, 'FontSize', 14, 'Location', 'best');
set(gca, 'FontSize', 14);
title(sprintf('lambda = %dnm, h = %d, Air',ida,h));
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the solution(Optional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;
figure;
model.result('pg1').set('data', indicator.dset);
model.result('pg1').set('solnum', 4);

model.result('pg1').feature('surf1').set('expr', "log(ewfd.normE)");  
model.result('pg1').feature('surf1').set('descr', '|E|');

model.result('pg1').feature('surf1').set('colorlegend', 'on');
model.result('pg1').feature('surf1').set('rangecoloractive', 'on');
model.result('pg1').feature('surf1').set('rangedatamin', '-5');
mphplot(model, 'pg1');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh Sweeping COnvergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = [10];
h = [2];
delta_list = [7,8,9,10,11,12,13,14,15];
mesh_list = [50,60,70,80,90,100];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize dictionary (Map)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = containers.Map();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
model.param.set('initial_guess', num2str(n_initial_guess, '%.5f'));
model.param.set('mesh_size', '80[nm]');
model.param.set('mesh_size2', '100[nm]');
model.study('std1').feature('mode').set('neigs', '4');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Record Computation Time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current = 1;
total = length(width_list)*length(mesh_list)*length(delta_list);
startTime = datetime('now');  % Record start time
disp(['Start time: ', datestr(startTime, 'yyyy-mm-dd HH:MM:SS')]);



% We want kappa(mesh)
temp = zeros(length(mesh_list),1);

j = 1;
for delta = delta_list
    key = sprintf('d_%d',delta);
    result(key) = 0 ;%initialize it
    for m = mesh_list
        disp(sprintf('%d / %d',current,total));
        % Load
        model.param.set('core_width', [num2str(w, '%.3f'),'[um]']);
        model.param.set('core2_width', [num2str(w, '%.3f'),'[um]']);
        model.param.set('core_height', [num2str(h, '%.3f'),'[um]']);
        model.param.set('delta', [num2str(delta, '%.3f'),'[um]']);
        model.param.set('mesh_size2', [num2str(m, '%.2f'),'[nm]']);
        model.param.set('mesh_size', [num2str(m, '%.2f'),'[nm]']);
    
    
        % Run
        indicator = flag(333); 
        studyName = indicator.std; 
        model.study(studyName).run();
        
        % Read
        beta = mphglobal(model, '(ewfd.beta)', 'dataset', indicator.dset);
        ex = mphglobal(model, '(ex)', 'dataset', indicator.dset);
        [~, idx] = sort(ex, 'descend');
        n1 = idx(1);
        n2 = idx(2);
        kappa = abs(beta(n1) - beta(n2));
    
    
        % Save
        temp(j,1) = kappa;
        j = j+1; 
        current = current+1;
    end
    result(key) = temp;
end

endTime = datetime('now');  % Record end time
disp(['End time: ', datestr(endTime, 'yyyy-mm-dd HH:MM:SS')]);
elapsedTime = endTime - startTime;
disp(['Elapsed time: ', char(elapsedTime)]);
%%
convergence = result;
%delta_list = [7,8,9,10,11,12,13,14,15]
delta = 10;
m = find(delta_list == delta, 1);
key = sprintf('d_%d',delta);
figure;
y = result(key);
y = y(1+(m-1)*6:m*6);
scatter(mesh_list,y);


%%




%% Saving
mphsave(model, savePath);
disp(['Model saved to: ', savePath]);
%% Exit 
% ModelUtil.disconnect();
ModelUtil.remove('model');
disp('Model removed from memory.');
%%
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

function loader(model,ida,n_initial_guess,core_width,core_height)
    model.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
    model.param.set('initial_guess', num2str(n_initial_guess, '%.5f'));
    model.param.set('core_width', [num2str(core_width, '%.4f'),'[um]']);
    model.param.set('core_height', [num2str(core_height, '%.4f'),'[um]']);
    model.param.set('mesh_size', 'min(ida/10,100[nm])');
    disp('Configured');
end