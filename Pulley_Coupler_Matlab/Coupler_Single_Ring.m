%% Initialization
addpath('C:\Program Files\COMSOL\COMSOL55\Multiphysics\mli')
import com.comsol.model.*
import com.comsol.model.util.*
format long;
mphstart(2036); 
%%
modelPath1 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\PG_Coupler_Jan_9_2025\Single Ring_Body.mph";
savePath1 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\PG_Coupler_Jan_9_2025\Single Ring_Body.mph";

modelPath2 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\PG_Coupler_Jan_9_2025\Single Ring_Bus.mph";
savePath2 = "C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\PG_Coupler_Jan_9_2025\Single Ring_Bus.mph";

model1 = ModelUtil.create('Model1');
model1 = mphload(modelPath1, 'Model1');

model2 = ModelUtil.create('Model2');
model2 = mphload(modelPath2,'Model2');  

ModelUtil.showProgress(true);

%% Navigator
mphnavigator(model2);
    % Std1 is for small Region ewfd2
    % Std2 is for large Region ewfd
%% Global Parameters
ida_list = [532,1064];
radius_list  = [1250,1300,1350,1400];
scale = length(ida_list)*length(radius_list);
jj = 1;

for radius = radius_list
for ida = ida_list
nClad = 1.450721; %Sub Thox
% 1.481681554; %2% Long; 1064
% 1.498419381; %2% Long; 532


if ida > 800
    nCore = 1.481681554 ;
    mesh_bus = 100;
    mesh_body = 100;
else
    nCore = 1.498419381 ;
    mesh_bus = 50;
    mesh_body = 75;
end

pml_extend = 4;
width_body = 9;
width_bus = 2;
height = 2; 
gap_distance = 1.5;
radius2 = radius+gap_distance+width_bus/2+width_body/2;


sim_width = 5;
% Configure Parameters for Body Ring
model1.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
model1.param.set('core_mesh_size', [num2str(mesh_body, '%.2f'),'[nm]']);
model1.param.set('nClad', num2str(nClad, '%.5f'));
model1.param.set('nCore', num2str(nCore, '%.5f'));
model1.param.set('pml_extend', [num2str(pml_extend, '%.2f'),'[um]']);
model1.param.set('radius', [num2str(radius, '%.3f'),'[um]']);
model1.param.set('width', [num2str(width_body, '%.4f'),'[um]']);
model1.param.set('height', [num2str(height, '%.4f'),'[um]']);
disp('Body Configured')
% Configure Parameters for Bus Ring
model2.param.set('ida', [num2str(ida, '%.3f'),'[nm]']);
model2.param.set('core_mesh_size', [num2str(mesh_bus, '%.2f'),'[nm]']);
model2.param.set('nClad', num2str(nClad, '%.3f'));
model2.param.set('nCore', num2str(nCore, '%.3f'));
model2.param.set('radius', [num2str(radius, '%.3f'),'[um]']);
model2.param.set('pml_extend', [num2str(pml_extend, '%.2f'),'[um]']);
model2.param.set('width', [num2str(width_bus, '%.4f'),'[um]']);
model2.param.set('height', [num2str(height, '%.4f'),'[um]']);
model2.param.set('interest_width', [num2str(sim_width, '%.4f'),'[um]']);
disp('Bus Configured')

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
    iTE = top2Indices(1);
    iTM = top2Indices(2);
else
    iTE = top2Indices(2);
    iTM = top2Indices(1);
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
figure;
model1.result('pg1').set('data', indicator.dset);
model1.result('pg1').set('solnum', num2str(solnum));
model1.result('pg1').feature('surf1').set('expr', indicator.normE); % Ensure 'surf1' exists in pg1
model1.result('pg1').feature('surf1').set('descr', '|E|');  
mphplot(model1, 'pg1');
drawnow;

if ida > 800
    %  Outer Loop for Different Bus Widths
    % Scan Range for Bus Width
    width_bus_range = [2,4]; 
    % Scan Range for Gap Distance
    gap_i = 4;
    gap_f = 6.6;
    step_size = 0.5;
    gap_range = gap_i:step_size:gap_f;
    num_points = length(gap_range);
else 
    %  Outer Loop for Different Bus Widths
    % Scan Range for Bus Width
    width_bus_range = [0.5,2]; 
    % Scan Range for Gap Distance
    gap_i = 1;
    gap_f = 4;
    step_size = 0.5;
    gap_range = gap_i:step_size:gap_f;
    num_points = length(gap_range);
end

results = struct(); % Store results for each bus width
all = length(gap_range)*length(width_bus_range);

for width_bus = width_bus_range
    disp(['Exploring bus width = ', num2str(width_bus), ' um']);

    % Initialize Parameters for Current Bus Width
    gap_distance = 1.5;
    radius2 = radius + gap_distance + width_bus / 2 + width_body / 2;

    % Preallocate Arrays
    bus_TE = zeros(1, num_points);
    bus_TM = zeros(1, num_points);

    % Inner Loop for Gap Distance
    i = 1;
    for gap_distance = gap_range
        disp(['     Running gap_distance = ', num2str(gap_distance), ' um ',num2str(jj/all*100/scale),'%']);
        
        % Update Parameters for Current Gap Distance
        radius2 = radius + gap_distance + width_bus / 2 + width_body / 2;
        sim_width = 1.5;

        % First Run: Eigenvalue Calculation
        model2.param.set('initial_guess', ['nCore*(radius+0[um])']);
        model2.param.set('interest_width', [num2str(sim_width, '%.4f'), '[um]']);
        model2.param.set('radius', [num2str(radius2, '%.4f'), '[um]']);

        % Run the study
        model2.study(studyName).run();

        % Evaluate Core-to-Total Energy Ratio
        [rate, ieig] = max(mphglobal(model2, 'core_energy2/total_energy2', 'dataset', indicator.dset));
        %disp(['     In-Core ratio is ', num2str(rate)]);
        if rate < 0.5
            % Skip invalid runs
            bus_TE(i) = -1;
            bus_TM(i) = -1;
            disp(['Failed at gap_distance = ', num2str(gap_distance),' Set neff = r*n']);
            model2.param.set('initial_guess', ['nCore*(radius+0[um])']);
        else
        % Extract Eigenvalue
        eigs = mphglobal(model2, 'real(ewfd2.neff)', 'dataset', indicator.dset);
        initial_guess = eigs(ieig);
        end
        % Second Run: Momentum Calculation
        sim_width = 5;
        model2.param.set('initial_guess', num2str(initial_guess, '%.6f'));
        model2.param.set('interest_width', [num2str(sim_width, '%.4f'), '[um]']);
        model2.study(studyName).run();

        % Extract TE and TM Values
        [~, iTE] = max(mphglobal(model2, 'ez2', 'dataset', indicator.dset));
        [~, iTM] = max(mphglobal(model2, 'ex2', 'dataset', indicator.dset));
        check = mphglobal(model2, 'core_energy2/total_energy2', 'dataset', indicator.dset);

        if check(iTE) > 0.7
            bus_TE(i) = mphglobal(model2, '(ewfd2.beta)', 'dataset', indicator.dset, 'solnum', iTE);
        else
            bus_TE(i) = -1;
        end

        if check(iTM) > 0.7
            bus_TM(i) = mphglobal(model2, '(ewfd2.beta)', 'dataset', indicator.dset, 'solnum', iTM);
        else
            bus_TM(i) = -1;
        end
        i = i + 1;
        jj = jj+1;
    end
    % Replace decimal point in width_bus with underscore for field name
    field_name = ['w', strrep(num2str(width_bus), '.', '_')];
    results.(field_name) = struct('gap_range', gap_range, 'bus_TE', bus_TE, 'bus_TM', bus_TM);
end
results.('body_TE') = body_TE;
results.('body_TM') = body_TM;
disp('Done exploring all bus widths.');
% Plot
figure;
for width_bus = width_bus_range
    % Get the results for the current bus width
    field_name = ['w', strrep(num2str(width_bus), '.', '_')];
    data = results.(field_name);
    disp(data);
    % Plot bus_TE for the current bus width
    plot(data.gap_range, data.bus_TE, '-o', 'DisplayName', ['bus\_TE (width = ', num2str(width_bus), ')'], 'LineWidth', 1.5);
    hold on;
end

% Plot the horizontal line for body_TE (black dashed line)
yline(body_TE, '--k', 'DisplayName', 'body\_TE', 'LineWidth', 1.5);

% Optionally, plot bus_TM lines if needed
if 0
    for width_bus = width_bus_range
        field_name = ['w', strrep(num2str(width_bus), '.', '_')];
        data = results.(field_name);
        plot(data.gap_range, data.bus_TM, '-x', 'DisplayName', ['bus\_TM (width = ', num2str(width_bus), ')'], 'LineWidth', 1.5);
        hold on;
    end
    % Plot the horizontal line for body_TM (green dashed line)
    yline(body_TM, '--g', 'DisplayName', 'body\_TM', 'LineWidth', 1.5);
end

% Labeling and grid setup
xlabel('Gap Distance', 'Interpreter', 'latex');
ylabel('$\beta_{eff}$','Interpreter', 'latex');
legend('show', 'Location', 'best');
grid on;
title(sprintf(['TE Modes vs. Gap Distance @%.1fnm \nRing %.1f' ...
    'um , %.2fx%.2f um'], ida,radius,width_body, height), 'Interpreter', 'latex');
hold off;

% Base folder path
base_folder_path = 'C:\Users\Dirk\Desktop\Hongrui_Yan_Simulation\PG_Coupler_Jan_9_2025\log';

% Create a folder name based on ring dimensions
folder_name = sprintf('Ring%.1fum_%.2fx%.2fum', radius, width_body, height);

% Combine the base path with the folder name
new_folder_path = fullfile(base_folder_path, folder_name);

% Ensure the ring dimension folder exists
if ~exist(new_folder_path, 'dir')
    mkdir(new_folder_path);
end

% File name for results
results_file_name = sprintf('@%.1fnm_Ring%.1fum_%.2fx%.2fum.mat', ida, radius, width_body, height);
results_file_path = fullfile(new_folder_path, results_file_name);

% Save the results structure as a MAT file
save(results_file_path, 'results');
disp(['Results saved as ', results_file_path]);

% File name for figure
figure_file_name = sprintf('@%.1fnm_Ring%.1fum_%.2fx%.2fum.fig', ida, radius, width_body, height);
figure_file_path = fullfile(new_folder_path, figure_file_name);

% Save the current figure as a FIG file
savefig(figure_file_path);
disp(['Figure saved as ', figure_file_path]);
end
end
%% Saving
if 1
    mphsave(model2, savePath2);
    disp(['Model saved to: ', savePath2]);
    mphsave(model1, savePath1);
    disp(['Model saved to: ', savePath1]);
    % Exit 
    % ModelUtil.disconnect();
    ModelUtil.remove('model');
    disp('Model removed from memory.');
end
%%
function output = flag(x)
   % x = 0, means small region, ewfd2
   output = struct();
   if x==0
       output.dset = 'dset1';
       output.normE = 'ewfd2.normE';
       output.std = 'std1';
   else
       output.dset = 'dset2';
       output.normE = 'ewfd1.normE';      
       output.std = 'std2';
   end 
end

function D_coeffs=dispersion(omega,beta,lambda_0)
    c=299792458;
    lambda_pump=lambda_0*1e-9;
    freq=omega;
    omega_pump=2*pi*c/lambda_pump;
    modenumber1=beta;%beta_eff is mu
    m01=interp1(freq,modenumber1,omega_pump); %find mu at pump frequency
    f2=polyfit(modenumber1-m01,omega-omega_pump,5);
    D1=f2(5);
    D2=2*f2(4);
    D3=6*f2(3);
    D4=24*f2(2);
    D5=120*f2(1);

    D_coeffs = [D1, D2, D3, D4, D5];

    D_int=omega-omega_pump-D1.*(modenumber1-m01);
    D_int_cal=1/2*D2*(modenumber1-m01).^2;
    figure;
    plot(modenumber1 - m01, D_int, 'x--', 'DisplayName', 'w - D_1  \mu', 'LineWidth', 1.5,'MarkerSize',8,'Color','black'); 
    hold on;
    plot(modenumber1-m01,D_int_cal,'DisplayName',  '1/2 D_2 \mu^2','LineWidth',1.5,'Color','blue')
    legend('Location', 'best');
    legend('show');
    FSR = D1 / (2 * pi)/(1e9); % FSR from D1
    b = D2 / (2 * pi); % FSR from D2
    title(sprintf('FSR = D_1 / 2\\pi = %.3f Ghz, D_2 / 2\\pi = %.3f', FSR, b));
end


%{
indicator = flag(0); %Do Small Region ewfd2
studyName = indicator.std; 


% Scan Range
x=8;
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
%% Extract Info
eigs = mphglobal(model, '(ewfd2.neff)', 'dataset', indicator.dset);
disp(eigs);
%%
dispersion(omega_mu,mu,1064)
%}