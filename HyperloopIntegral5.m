close all

% ALL NUMBERS ARE IN METRIC AND ARE CONVERTED TO FEET WHERE NECESSARY

% Time (Defined arbitrarily to get numbers but has no effect on
%       calculation)
t = 1;


% Physical Constants
rho = 3.99e-8;                  % Track resistivity
g_c = 0.413 * 0.0254;           % Track thickness
rho_s = rho / g_c;              % Track surface resistivity
mu_0 = pi * 4e-7;               % Vacuum permeability
rho_c = 1.72e-8;                % Coil resisitivity
sigma = 5.67e-8;                % Boltzmann's constant


% Stability's Constants
h_nz1 = 0.7 * 0.0254;
h_nz2 = (5-.413-.7) * 0.0254;
h_nz3 = 0.5 * 0.0254;
sigma_sh = 0.5 * 0.0254;        % Minimum available horizontal variation
w_nz1 = 0.75;
w_nz2 = 0.412;
sigma_sw = 0.5 * 0.0254;        % Minimum available vertical variation


%########### Straddle LIM Configurations #############%
%{
%--------------Config 1 (60N) (Closest)---------------%
g = g_c + (2 * (sigma_sh + 0)); 
W = (3.42 - w_nz1 - (2*sigma_sw))* 0.0254;
%}
%{
%--------------Config 2 (40N)-------------------------%
%g = g_c + (2 * (sigma_sh + h_nz1)); 
%W = (3.42 - (2*sigma_sw))* 0.0254;
%}
%{
%--------------Config 3 (6N)--------------------------%
%g = g_c + (2 * (sigma_h + h_nz1 + h_nz2));
%W = 3.42 * 0.0254;
%}
%{
%--------------Config 4 (6N)--------------------------%
%g = g_c + (2 * (sigma_sh + h_nz1 + h_nz2 + h_nz3)); 
%W = (3.42 + w_nz2) * 0.0254;
%}
%#####################################################%


%############## Top LIM Configurations ###############%
%--------------Config 5 (60N) (Closest)---------------%
g = g_c + (2 * (sigma_sh)); 
W = 5 * 0.0254;

%{
%--------------Config 6 (40N)-------------------------%
%g = g_c + (2 * (1.5*sigma_sh)); 
%W = 5 * 0.0254;
%}
%#####################################################%


% LIM Configuration Constants
L = 55 * 0.0254;                % LIM length
f = 10;                         % Starting frequency
r_w = 7.12e-4;                  % Wire radius
r_c = 14*r_w;                   % Coil radius
N = 127;                        % Windings per coil
T_w = 155 + 273;                % Wire temperature threshold
p = 4;                          % Number of poles (not pole pairs)
tau = L/(4*p);                  % Pole pitch

% Check for configuration validity
if r_c > (tau/3)
   disp("Coil width larger than tooth");
   return
end


% Desired LIM Qualities
max_speed = 111.76;
max_freq = (max_speed/(2*tau));


% Stuff to make the plot look nice
max_thrust_plot = 0;
max_freq_plot = 0;
max_thrust_speed = 0;
optimal_thrust = [];
optimal_freq = [];


% Resolution 
freq_increment = 1;             % Frequency step size 
v_fidelity = 100;               % Resolution of each graph


% This v_s and the second condition in the while loop are to ensure no
% runaway simulations
v_s = 0;
while (max_thrust_speed < max_speed) && (v_s < max_speed + 100)
    
    clf
    
    % Other Parameters
    v_s = 2*f*tau;                              % Motor speed
    v = linspace(1,v_s, v_fidelity);            % Velocity array
    omega = 2*pi*f;                             % Angular frequency
    c = (4*rho_s*g*omega)./(mu_0*(v.^2));       % Intermediate term
    tau_e = (v*pi)/(omega*(2^0.5))...           % Shockwave half wavelength
        .*(1 + (1 + (c.^2)).^0.5).^0.5;
    beta = (rho_s*g)./(v*mu_0);                 % Intermediate term
    I = (2*pi^2*r_w^3*sigma*T_w^4/rho_c)^0.5;   % Max tolerable current
    J_1 = (3*(2^0.5)*N*I)/(p*tau);              % LIM surface current
    
    
    
    %######## Magnetic Field Amplitudes ########%
    
    % LIM's Magnetic Field
    B_s = J_1 ./ (((pi*g)/(tau*mu_0))^2 ...     % LIM magnetic field magnitude
        + ((1/rho_s)*(v_s-v)).^2).^0.5;
    delta_s = atan((pi*rho_s*g)...              % LIM magnetic field phase
        ./(mu_0*tau*(v_s-v)));
    
    % Intermediate terms for B1 and B2
    a = ((mu_0*v)/(rho_s*g)).^2;
    b = (4*omega*mu_0)/(rho_s*g);
    X = ((((a.^2 + b^2).^0.5) + a)./2).^0.5;
    
    
    % Front Magnetic Shockwave
    alpha_1 = (2*rho_s*g)...                    % Front magnetic shockwave decay constant
        ./((rho_s*g.*X) - (mu_0.*v));
    
    % Front magnetic shockwave intermediate terms
    gamma_1 = ((1 + (beta/alpha_1)).^2 + ((beta*pi)/tau_e).^2).^(-1);
    term1 = -(rho_s.*J_1./v) - B_s.*cos(delta_s) + (beta.*pi.*B_s/tau).*sin(delta_s);
    term2 = 1 + (beta/alpha_1);
    term3 = (beta.*pi.*B_s/tau).*cos(delta_s) + B_s.*sin(delta_s);
    term4 = beta*pi/tau;
    
    delta_1 = atan((term1.*term4 ...            % Front magnetic shockwave phase
        + term3.*term2)/(term1.*term2...
        + term3.*term4));
    B_1 = (term1.*term4 + term3.*term2)...      % Front magnetic shockwave amplitude
        *gamma_1/cos(delta_1);
    
    
    % Rear Magnetic Shockwave
    alpha_2 = (2*rho_s*g)...                    % Rear magnetic shockwave decay constant
        ./((rho_s*g.*X) + (mu_0.*v));
    
    % Rear magnetic shockwave intermediate terms
    gamma_2 = ((1 - (beta/alpha_2)).^2 + ((beta*pi)./tau_e).^2).^(-1);
    term5 = B_1*cos(delta_1);
    term6 = 1 - (exp(-L./alpha_1).*cos(pi.*L./tau_e));
    term7 = B_1*sin(delta_1);
    term8 = exp(-L./alpha_1).*sin(pi.*L./tau_e);
    term9 = -1 - (beta/alpha_1);
    term10 = 1 + (beta/alpha_2);
    term11 = -(beta*pi)/tau_e;
    term12 = (beta*pi)/tau_e;
    term13 = term5.*term6 - term7.*term8;
    term14 = term5.*term8 + term7.*term6;
    term15 = term9.*term10 - term11.*term12;
    term16 = term9.*term12 + term10.*term11;
    term17 = gamma_2.*(term13.*term15 - term14.*term16);
    term18 = gamma_2.*(term15.*term14 + term13.*term16);
    
    delta_2 = atan(term18./term17);             % Rear magnetic shockwave phase
    B_2 = gamma_2.*term17./cos(delta_2);        % Rear magnetic shockwave amplitude
    
    %###########################################%
    
    
    % Total Magnetic Field
    func = @(x) J_1.*B_s.*cos(delta_s) ...
        + J_1.*B_1.*exp(-x./alpha_1).*cos((pi.*x).*((1./tau)-(1./tau_e))+delta_1) ...
        + J_1.*B_2.*exp(-x./alpha_2).*cos((pi.*x).*((1./tau)+(1./tau_e))+delta_2);
    
    % Total Thrust
    T = W*integral(func,0,L,'ArrayValued',true);
    
    
     Potential for frequency vs speed plot
    {
    tiledlayout(2, 1);
    set(gcf, 'Position',  [100, 100, 1000, 800])
    nexttile
    }
    
    % Makes the plots look nice
    max_T = find(T == max(T(:)));
    if v(max_T) > max_thrust_speed
       max_thrust_speed = v(max_T); 
    end
    if max(T)>max_thrust_plot
       max_thrust_plot = max(T); 
    end
    
    % Optimal Thrust vs Speed Curve
    optimal_thrust = [optimal_thrust; max_thrust_speed max(T)];
    
    
    %Plots
    hold on
    plot(v, T);
    xline(v_s, '--');
    plot([max_thrust_speed, max_thrust_speed], [0, max(T)]);
    plot(optimal_thrust(:,1), optimal_thrust(:,2));
    xlabel('Speed (mph)');
    ylabel('Thrust (Newtons)');
    axis([0 v_s 0 max_thrust_plot]);
    hold off
    
    
    % Potential for frequency vs speed plot
    %{
    %nexttile
    %optimal_freq = [optimal_freq; max_thrust_speed f];
    %if max(optimal_freq(:,2))>max_freq_plot
       %max_freq_plot = max(optimal_freq(:,2)); 
    %end
    
    plot(optimal_freq(:,1), optimal_freq(:,2));
    xlabel('Speed (mph)');
    ylabel('Frequency (Hz)');
    axis([0 v_s 0 max_freq_plot]);
    %}
    
    
    f = f + freq_increment;                     % Incrementing frequency
    
    % Comment this out if you don't want to draw the whole simulation
    drawnow    
end