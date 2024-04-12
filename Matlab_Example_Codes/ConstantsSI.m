%% Constants in SI units

SI = struct();

SI.h = 6.62607004e-34;              % Planck's constant / J.s
SI.hbar = SI.h/(2*pi);              % Reduced Planck's constant / J.s
SI.kb = 1.38064852e-23;             % Boltzmann's constant / J.K^(-1)
SI.Na = 6.022140857e23;             % Avogadro number / mol-1
SI.c = 2.997925e8;                  % Speed of light / m.s-1
SI.mu0 = 4*pi * 1e-7;               % Vacuum permeability / N.A^(-2)
SI.eps0 = 1/(SI.mu0*SI.c^2);        % Vacuum permittivity / J^(-1).C^2.m^(-1)
SI.e = 1.6021766208e-19;            % Elementary charge / C
SI.Faraday = SI.Na*SI.e;            % Faraday constant / C.mol^(-1)
SI.m_e = 9.10938e-31;               % Electron mass / kg
SI.m_n = 1.674927351e-27;           % Neutron mass / kg
SI.m_p = 1.67262171e-27;            % Proton mass / kg
SI.m_u = 1e-3/SI.Na;                % Atomic mass unit / kg
SI.R = SI.Na*SI.kb;                 % Ideal gas constant / J.K^(-1).mol^(-1)
SI.muB = 0.5*SI.e*SI.hbar/SI.m_e;   % Bohr magneton / J.T^-1
SI.muN = 0.5*SI.e*SI.hbar/SI.m_p;   % Nuclear magneton / J.T^-1

SI.a0 = 5.291772106712e-11;       % Length unit (Bohr) / m
SI.Eh = 4.35974465054e-18;        % Energy unit (Hartree) / J
SI.ta = 2.41888432650516e-17;     % Time unit (hbar_/Hartree) / s
SI.Ta = 3.157746455e5;            % Temperature unit (Hartree/kb_) / K



