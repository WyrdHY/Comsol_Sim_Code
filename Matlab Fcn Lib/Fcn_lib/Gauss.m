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
