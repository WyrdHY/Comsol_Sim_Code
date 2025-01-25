%% Initialization
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