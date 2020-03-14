#!/usr/bin/env python3
# ALL NUMBERS ARE IN METRIC AND ARE CONVERTED TO FEET WHERE NECESSARY

# Time (Defined arbitrarily to get numbers but has no effect on calculation)

# This file was converted from MATLAB to Python by Controls

import math
import argparse
import matplotlib.pyplot as plt
import numpy as np

class Propulsion:
    # Default parameters are physical and stability's constants
    def __init__(
        self, 
        h_nz1 = 0.7 * 0.0254,
        h_nz2 = (5-.413-.7) * 0.0254,
        h_nz3 = 0.5 * 0.0254,
        sigma_sh = 0.5 * 0.0254,
        w_nz1 = 0.75,
        w_nz2 = 0.412,
        sigma_sw = 0.5 * 0.0254,
        rho = 3.99e-8,
        g_c = 0.413 * 0.0254,
        rho_c = 1.72e-8):

        # Universal Physical Constants
        self.mu_0 = math.pi * 4e-7          # Vacuum permeability
        self.sigma = 5.67e-8                # Boltzmann's constant

        # Stability's Constants
        self.h_nz1 = h_nz1
        self.h_nz2 = h_nz2
        self.h_nz3 = h_nz3
        self.sigma_sh = sigma_sh    # Minimum available horizontal variation
        self.w_nz1 = w_nz1
        self.w_nz2 = w_nz2
        self.sigma_sw = sigma_sw    # Minimum available verical variation  

        # Physical Constants
        self.rho = rho                      # Track resistivity
        self.g_c = g_c                      # Track thickness
        self.rho_s = self.rho / self.g_c    # Track surface resistivity
        self.rho_c = rho_c                  # Coil resisitivity
        
        # Top LIM Config 5 (60N) (Closest)
        self.g = g_c + (2 * (sigma_sh))
        self.W = 5 * 0.0254

        # LIM Configuration Constants
        self.L = 55 * 0.0254                 # LIM length
        self.f = 10                          # Starting frequency
        self.r_w = 7.12e-4                   # Wire radius
        self.r_c = 14*self.r_w               # Coil radius
        self.N = 127                         # Windings per coil
        self.T_w = 155 + 273                 # Wire temperature threshold
        self.p = 4                           # Number of poles (not pole pairs)
        self.tau = self.L/(4*self.p)         # Pole pitch

    def check_configuration_validity(self):
        if self.r_c > (self.tau/3):
            print("Coil width larger than tooth")

    def setup_simulation(self, max_speed = 111.76):
        # Desired LIM Qualities
        self.max_speed = max_speed
        self.max_freq = (self.max_speed/(2*self.tau))

        # Stuff to make the plot look nice
        self.max_thrust_plot = 0
        self.max_freq_plot = 0
        self.max_thrust_speed = 0
        self.optimal_thrust = []
        self.optimal_freq = []

        # Resolution
        self.freq_increment = 1
        self.v_fidelity = 100

def simulate():
    # v_s and the second condition in the while loop prevent runaway simulations
    v_s = 0
    f = 10
    max_thrust_speed = 0
    max_speed = 111.76
    freq_increment = 1
    while (max_thrust_speed < max_speed and v_s < max_speed + 100):
        # Other Parameters
        max_thrust_speed += 1
        print(v_s)
        print(f)
        v_s = 2 * f * tau                           # Motor speed             
        v = np.linspace(1, v_s, num=v_fidelity)     # Velocity array
        omega = 2 * math.pi * f                     # Angular frequency
        c = (4*rho_s*g*omega) / (mu_0*(v**2))       # Intermediate term
        tau_e = (v*math.pi)/(omega*(2**0.5)) \
        * (1 + (1 + (c**2))**0.5)**0.5             # Shockwave half wavelength
        beta = (rho_s*g) / (v*mu_0)                 # Intermediate term
        I = (2*math.pi**2*r_w**3*sigma*T_w \
        **4/rho_c)**0.5                             # Max tolerable current
        J_1 = (3*(2**0.5)*N*I)/(p*tau)              # LIM surface current
    
    
        # %######## Magnetic Field Amplitudes ########%
        
        # % LIM's Magnetic Field
        B_s = J_1 / (((math.pi*g)/(tau*mu_0))**2 /  # LIM magnetic field magnitude
            + ((1/rho_s)*(v_s-v))**2)**0.5
        delta_s = np.arctan((math.pi*rho_s*g) \
            / (mu_0*tau*(v_s-v)))                   # LIM magnetic field phase
        
        # % Intermediate terms for B1 and B2
        a = ((mu_0*v)/(rho_s*g))**2
        b = (4*omega*mu_0)/(rho_s*g)
        X = ((((a**2 + b**2)**0.5) + a) / 2)**0.5
        
        
        # % Front Magnetic Shockwave
        alpha_1 = (2*rho_s*g) \
            / ((rho_s * g * X) - (mu_0 * v))   # Front magnetic shockwave decay constant
        
        # % Front magnetic shockwave intermediate terms
        
        gamma_1 = ((1 + (beta/alpha_1))**2 + ((beta*math.pi)/tau_e)**2)**(-1)
        #TODO: fix term1 and term3
        term1 = -(rho_s * np.divide(J_1, v, np.vectorize(np.int))) - B_s * math.cos(delta_s) + (beta * math.pi *B_s/tau) * math.sin(delta_s)
        term2 = 1 + (beta/alpha_1)
        # term3 = (beta * math.pi * B_s / tau) * math.cos(delta_s) + B_s * math.sin(delta_s)
        term4 = beta*math.pi/tau
        # delta_1 = np.arctan((term1 * term4 \
        #     + term3 * term2)/(term1 * term2 \
        #     + term3 * term4))             # Front magnetic shockwave phase
        # B_1 = (term1 * term4 + term3 * term2) \
        #     * gamma_1/math.cos(delta_1)         # Front magnetic shockwave amplitude
        
        
        # Rear Magnetic Shockwave
        alpha_2 = (2 * rho_s * g) \
            / ((rho_s * g * X) + (mu_0 * v))        # Rear magnetic shockwave decay constant
        
        # Rear magnetic shockwave intermediate terms
        gamma_2 = ((1 - (beta/alpha_2))**2 + ((beta*math.pi) / tau_e)**2)**(-1)
        term5 = B_1 * math.cos(delta_1)
        term6 = 1 - (exp(-L / alpha_1) * math.cos(math.pi * L / tau_e))
        # term7 = B_1*sin(delta_1);
        # term8 = exp(-L./alpha_1).*sin(pi.*L./tau_e);
        # term9 = -1 - (beta/alpha_1);
        # term10 = 1 + (beta/alpha_2);
        # term11 = -(beta*pi)/tau_e;
        # term12 = (beta*pi)/tau_e;
        # term13 = term5.*term6 - term7.*term8;
        # term14 = term5.*term8 + term7.*term6;
        # term15 = term9.*term10 - term11.*term12;
        # term16 = term9.*term12 + term10.*term11;
        # term17 = gamma_2.*(term13.*term15 - term14.*term16);
        # term18 = gamma_2.*(term15.*term14 + term13.*term16);
        
        # delta_2 = atan(term18./term17);             % Rear magnetic shockwave phase
        # B_2 = gamma_2.*term17./cos(delta_2);        % Rear magnetic shockwave amplitude
        
        # %###########################################%
        
        
        # % Total Magnetic Field
        # func = @(x) J_1.*B_s.*cos(delta_s) ...
        #     + J_1.*B_1.*exp(-x./alpha_1).*cos((pi.*x).*((1./tau)-(1./tau_e))+delta_1) ...
        #     + J_1.*B_2.*exp(-x./alpha_2).*cos((pi.*x).*((1./tau)+(1./tau_e))+delta_2);
        
        # Total Thrust
        # T = W*integral(func,0,L,'ArrayValued',true)
        
        
        #  Potential for frequency vs speed plot
        # {
        # tiledlayout(2, 1);
        # set(gcf, 'Position',  [100, 100, 1000, 800])
        # nexttile
        # }
        
        # Makes the plots look nice
        # max_T = np.amax(T)
        # if v(max_T) > max_thrust_speed:
        #    max_thrust_speed = v(max_T)
        # if np.amax(T) > max_thrust_plot:
        #    max_thrust_plot = np.amax(T)
        
        # % Optimal Thrust vs Speed Curve
        # optimal_thrust = [optimal_thrust; max_thrust_speed max(T)];
        
        
        # %Plots
        # hold on
        # plot(v, T);
        # xline(v_s, '--');
        # plot([max_thrust_speed, max_thrust_speed], [0, max(T)]);
        # plot(optimal_thrust(:,1), optimal_thrust(:,2));
        # xlabel('Speed (mph)');
        # ylabel('Thrust (Newtons)');
        # axis([0 v_s 0 max_thrust_plot]);
        # hold off
        
        
        # % Potential for frequency vs speed plot
        # %{
        # %nexttile
        # %optimal_freq = [optimal_freq; max_thrust_speed f];
        # %if max(optimal_freq(:,2))>max_freq_plot
        #    %max_freq_plot = max(optimal_freq(:,2)); 
        # %end
        
        # plot(optimal_freq(:,1), optimal_freq(:,2));
        # xlabel('Speed (mph)');
        # ylabel('Frequency (Hz)');
        # axis([0 v_s 0 max_freq_plot]);
        # %}
        f = f + freq_increment                     # Incrementing frequency        

