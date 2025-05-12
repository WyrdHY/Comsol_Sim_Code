function output = flag(x)
   % x = 0, means small region, ewfd2
   % x = 333 means no need to distinguish them
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

   if x == 333
        output.dset = 'dset1';
        output.normE = 'ewfd1.normE';      
        output.std = 'std1';
   end
end
function n_interp = load_and_interpolate_n(filepath)
    % Load the data from the provided file path
    data = load(filepath);
    wavelength_nm = data(:,1);  % First column: wavelength (nm)
    if wavelength_nm(1) <0.1
        wavelength_nm = wavelength_nm*1e9
    end
    n_values = data(:,2);       % Second column: refractive index (n)

    % Create the interpolation function
    n_interp = @(wavelength) interp1(wavelength_nm, n_values, wavelength, 'spline');

    % Display a confirmation message
    disp('Interpolation function created successfully.');
    disp(['Wavelength range: ', num2str(min(wavelength_nm)), ' nm to ', num2str(max(wavelength_nm)), ' nm']);
end

%% Gaussian Mode 
function [E_gauss_x,E_gauss_y]=Gauss(w0,X,Y,polarization)
    % polarization = 1 means TE, Eup
    % polarization = 0 means TM, Eright
    x0 = 0; y0 = 0; A  = 1;  
    if ~polarization
        E_gauss_y = A * exp(-((X - x0).^2 + (Y - y0).^2)/(w0^2));
        E_gauss_x = zeros(size(E_gauss_y));
    else
        E_gauss_x = A * exp(-((X - x0).^2 + (Y - y0).^2)/(w0^2));
        E_gauss_y = zeros(size(E_gauss_x));
    end
end

function [E_gauss_x,E_gauss_y]=Ellptical_Gauss_z(X,Y,polarization,z,theta_x,theta_y,ida)
    % polarization = 1 means TE, Eup
    % polarization = 0 means TM, Eright
    % ida is in nm
    % theta is the divergence angle in deg
    % z is in um
    z = z*1e-6;
    ida = ida*1e-9;
    theta_x=theta_x*pi/180;
    theta_y=theta_y*pi/180;
    wx0 = ida/(theta_x*pi);
    wy0 = ida/(theta_y*pi);
    
    strechx = ida*z/(pi*wx0^2);
    strechy = ida*z/(pi*wy0^2);
    wxz = wx0 * sqrt(1+(strechx)^2);
    wyz = wy0 * sqrt(1+(strechy)^2);

    x0 = 0; y0 = 0; A  = wx0*wy0/(wxz*wyz);  

    E_field = A * exp(-((X - x0).^2 / wxz^2 + (Y - y0).^2 / wyz^2));

    if ~polarization
        E_gauss_y = E_field;
        E_gauss_x = zeros(size(E_gauss_y));
    else
        E_gauss_x = E_field;
        E_gauss_y = zeros(size(E_gauss_x));
    end
end





function compare_plot(w0, solnum, Gauss_polar, model, indicator)
    % Function to compare Gaussian beam with COMSOL field
    % w0: Beam waist in um
    % solnum: Solution number (TE or TM mode)
    % Gauss_polar: Polarization direction (0 for x, 1 for y)
    % model: COMSOL model object
    % indicator: Dataset structure
    
    % Parameters for Gaussian Beam
    nx = 500;  
    ny = 300;  
    w0 = w0*1e-6;
    % Meshing
    rx = max(2*w0, 10e-6);
    ry = max(2*w0, 10e-6);
    xVec = linspace(-rx, rx, nx);
    yVec = linspace(-ry, ry, ny);
    [X, Y] = meshgrid(xVec, yVec); 
    coords = [X(:)'; Y(:)'];
    
    % Extract COMSOL fields
    Ex_data = reshape(mphinterp(model, 'ewfd.Ex', 'coord', coords, 'dataset', indicator.dset, 'solnum', solnum), ny, nx);
    Ey_data = reshape(mphinterp(model, 'ewfd.Ey', 'coord', coords, 'dataset', indicator.dset, 'solnum', solnum), ny, nx);

    % Generate Gaussian beam field
    [E_gauss_x, E_gauss_y] = Gauss(w0, X, Y, Gauss_polar); 

    % Convert x and y to micrometers
    xVec_um = xVec * 1e6;
    yVec_um = yVec * 1e6;

    % Create figure with side-by-side subplots
    figure;
    colormap jet;  

    % Plot Gaussian field
    % Set dynamic title based on polarization
    if Gauss_polar == 1
        title_str = sprintf('w0 = %.2f \x03BCm, TE Mode', w0 * 1e6);
    else 
        title_str = sprintf('w0 = %.2f \x03BCm, TM Mode', w0 * 1e6);
    end
    subplot(1,2,1);
    imagesc(xVec_um, yVec_um, sqrt(E_gauss_y.^2+E_gauss_x.^2));
    
    title(title_str);
    set(gca, 'YDir', 'normal');  
    axis image;                 
    colorbar;
    xlabel('x [\mum]', 'FontSize', 12);
    ylabel('y [\mum]', 'FontSize', 12);

    % Plot COMSOL field
    subplot(1,2,2);
    if Gauss_polar == 1
        imagesc(xVec_um, yVec_um, abs(Ey_data));
    else 
        imagesc(xVec_um, yVec_um, abs(Ex_data));
    end
    set(gca, 'YDir', 'normal');  
    title('COMSOL Simulation')
    axis image;                  
    colorbar;




    xlabel('x [\mum]', 'FontSize', 12);
    ylabel('y [\mum]', 'FontSize', 12);

    % Overall figure title
    sgtitle("Compare")

end
%% Dispersion
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

%%
a =1;
disp(~a);