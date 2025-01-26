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
