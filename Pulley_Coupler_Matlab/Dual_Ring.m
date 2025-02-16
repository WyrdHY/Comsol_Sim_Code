addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Matlab Fcn Lib')));%% Initialization
addpath('C:\Program Files\COMSOL\COMSOL55\Multiphysics\mli')
import com.comsol.model.*
import com.comsol.model.util.*
format long;
mphstart(2036); 
%%
modelPath1 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\PG_Coupler_Jan_9_2025\Dual Ring.mph";
savePath1 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Model\PG_Coupler_Jan_9_2025\Dual Ring.mph";

model1 = ModelUtil.create('Model1');
model1 = mphload(modelPath1, 'Model1');

ModelUtil.showProgress(true);
%% Navigator
mphnavigator(model1);
    % Std1 is for small Region ewfd2
    % Std2 is for large Region ewfd
%% Global Parameters
ida=1550;
width = 30;
width2 = 5.93173617040035;
gap_distance=0.8;
radius=3296.747;
radius2 = radius + gap_distance + (width+width2)/2 ;
height = 4;
cladding_height = 5;
mesh_size = 100; %nm


n=load_and_interpolate_n("C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\Code\Material\Sellmeier Fitting\Sellmeier_2%_Long Anneal.txt");
nCore = n(ida);
%%
model1.param.set('nCore', [num2str(nCore, '%.5f')]);
model1.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
model1.param.set('radius', [num2str(radius, '%.3f'),'[um]']);
model1.param.set('radius2', [num2str(radius2, '%.8f'),'[um]']);
model1.param.set('width', [num2str(width, '%.3f'),'[um]']);
model1.param.set('width2', [num2str(width2, '%.8f'),'[um]']);
model1.param.set('cladding_height', [num2str(cladding_height, '%.3f'),'[um]']);
model1.param.set('height', [num2str(height, '%.3f'),'[um]']);
model1.param.set('core_mesh_size',  [num2str(mesh_size, '%.3f'),'[nm]']);
model1.param.set('initial_guess', 'nCore*radius2');

disp(['nCore ', num2str(nCore, '%.5f')]);
disp(['ida ', num2str(ida, '%.3f'), '[nm]']);
disp(['radius ', num2str(radius, '%.3f'), '[um]']);
disp(['radius2 ', num2str(radius2, '%.8f'), '[um]']);
disp(['width ', num2str(radius, '%.3f'), '[um]']);
disp(['width2 ', num2str(width2, '%.8f'), '[um]']);
disp(['cladding_height ', num2str(cladding_height, '%.3f'), '[um]']);
disp(['height ', num2str(height, '%.3f'), '[um]']);

disp('Successfully Configured');



%% Body Momentum
% Prepare Macro for Simulation 
indicator = flag(0); % Small Region ewfd2
studyName = indicator.std; 
% Body Momentum
disp('Running Body Momentum')
model1.study(studyName).run();
neff = mphglobal(model1, 'real(ewfd2.neff)', 'dataset', indicator.dset);
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
check = mphglobal(model1, 'core_energy2/total_energy2', 'dataset', indicator.dset);
if check(iTE)>0.7
    body_TE = mphglobal(model1, '(ewfd2.beta)', 'dataset', indicator.dset,'solnum',iTE);
else
    body_TE = -1;
end
if check(iTM)>0.7
    body_TM = mphglobal(model1, '(ewfd2.beta)', 'dataset', indicator.dset,'solnum',iTM);
else
    body_TM = -1;
end
disp('Body Momentum Done')
solnum = iTE;

%Updates the Initial Guess
initial_guess = neff(iTE);
model1.param.set('initial_guess', [num2str(initial_guess, '%.7f')]);

%% Large Region
indicator = flag(1); % Small Region ewfd2
studyName = indicator.std; 
% Body Momentum
disp('Running Body Momentum')
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
if check(iTE)>0.7
    body_TE = mphglobal(model1, '(ewfd.beta)', 'dataset', indicator.dset,'solnum',iTE);
else
    body_TE = -1;
end
if check(iTM)>0.7
    body_TM = mphglobal(model1, '(ewfd.beta)', 'dataset', indicator.dset,'solnum',iTM);
else
    body_TM = -1;
end
disp('Body Momentum Done')
solnum = iTE;
%% Visualize it
solnum =iTE;
figure;
model1.result('pg1').set('data', indicator.dset);
model1.result('pg1').set('solnum', num2str(solnum));
model1.result('pg1').feature('surf1').set('expr', '(ewfd.normE)'); % Ensure 'surf1' exists in pg1
model1.result('pg1').feature('surf1').set('descr', '|E|');  
model1.result('pg1').feature('arws1').set('expr', {'abs(ewfd.Er)','abs(ewfd.Ez)'});
mphplot(model1, 'pg1');
mphplot(model1, 'pg1', 'rangenum', 1);
drawnow;

%% Saving
if 1

    mphsave(model1, savePath1);
    disp(['Model saved to: ', savePath1]);
    % Exit 
    % ModelUtil.disconnect();
    ModelUtil.remove('model');
    disp('Model removed from memory.');
end




