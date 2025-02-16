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
