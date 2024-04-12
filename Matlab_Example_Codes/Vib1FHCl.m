%Reaction: HCl + F -> HF + Cl
clear;


% *** Propagation time parameters
%(in fs)
tstep  = 0.01;     % propagation step
tout   = 30;        % save wf every tout fs
tfinal = 1200;       % final propagation time

n_start = 0;
n_end   = 0;

n_steps = floor(tfinal/tstep);      % Total propagation steps

%% Set LaTeX interpretar as default
close all;
Conv = Converters;



%Colours
% Define RGB values (as fractions between 0 and 1)
Grey = [69, 69, 69]/255; % This represents a shade of grey
NavyBlue = [0, 0, 128]/255; % This represents a shade of NavyBlue
Bordeaux = [95, 2, 31]/255; % This represents a shade of Bordeaux

% LaTeX interpreter for text in plots
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultColorbarTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontName','Verdana');
set(groot,'defaultAxesFontSize', 18);


%% Initialization
ConstantsSI
Parameters

%Convert to atomic units
alphaHF = Conv.Froma0toA(alphaHF);
ReqHF   = Conv.FromAtoa0(ReqHF);
DeHF    = Conv.FromeVtoHartree(DeHF);
alphaFCl = Conv.Froma0toA(alphaFCl);
ReqFCl   = Conv.FromAtoa0(ReqFCl);
DeFCl    = Conv.FromeVtoHartree(DeFCl);
alphaHCl = Conv.Froma0toA(alphaHCl);
ReqHCl   = Conv.FromAtoa0(ReqHCl);
DeHCl    = Conv.FromeVtoHartree(DeHCl);

% *** Masses and vib. frequencies
m = [mF mH mCl]/SI.m_e;   % In a.u.
mu = [muHF,muHCl,muFCl]/SI.m_e;
mu2 = [m(1)*(m(2)+m(3))/(m(1)+m(2)+m(3)), m(2)*m(3)/(m(2)+m(3))];
%w = [4401 4401 4401]*4.55633E-6; % period = 10fs
w = [100 4401 4401]*4.55633E-6; % period = 10fs


% *** Grid
Rmin = [0.2   0.2]*1E-10/SI.a0  ;  % in au
Rmax = [70   3]   *1E-10/SI.a0  ;  % in au
N(1) = floor(2^(10))     ;  % no. of grid points
N(2) = ceil(N(1) *(Rmax(2)-Rmin(2) )/(Rmax(1) - Rmin(1) ) );


% Create grid
r1 = linspace(Rmin(1), Rmax(1), N(1));
r2 = linspace(Rmin(2), Rmax(2), N(2));
LR = Rmax-Rmin;
[R1, R2] = ndgrid(r1,r2);
VContours = zeros(size(R1));
p1 = ifftshift(-floor(N(1)/2):ceil(N(1)/2)-1) * 2*pi/LR(1);   % k-space
p2 = ifftshift(-floor(N(2)/2):ceil(N(2)/2)-1) * 2*pi/LR(2);
[P1, P2] = ndgrid(p1, p2);
LP = [max(p1)-min(p1), max(p2)-min(p2)];

% *** Initial state: Gaussian wavepacket
x0 = [42*1E-10/SI.a0, ReqHCl];

n0      = 1;


D_e     = DeHCl;                               % Depth of the potential well
beta    = alphaHCl;                               % Morse potential parameter
hbar    = 1;                                    % Planck's reduced constant
lambda  = sqrt(8*mu(2)*D_e)/(beta*hbar);           % Define lambda value 
k       = beta^(2)*2*D_e;
omega   = sqrt(k/mu(2));
y       =  lambda*exp(-beta.*(R2-x0(2)));

E_Vib = @(mm) beta*hbar*sqrt(2*D_e/mu2(2))*(mm+1/2) - ...
        beta^(2)*(hbar)^(2)/(2*mu2(2))*(mm+1/2)^(2);
E_Vib0 = E_Vib(n0);


% Number of steps
n_silentSteps = round(tout/tstep);  % Write every n_silentSteps
n_outSteps = fix(n_steps/n_silentSteps);    % Number of output


%% Initialisation of the wavefunction

alpha    = (lambda - 2 * n0 - 1)/2;
Norm     = sqrt(beta * 2 * alpha*gamma(n0+1)/ ...
           (gamma(lambda - n0)));
Laguerre = laguerreL(n0, 2*alpha, y);


Wavefunction = Norm.* exp(-y/2).* Laguerre .* y.^(alpha);


sigma0 = [14.5337,    0.1654,    0.1654];

existingData = readmatrix('Vib1FHCl.txt');

for iii = n_start:n_end
%flux definition
flux = zeros(N(1), n_steps);
FluxIndex2 = ceil(N(2)*10/20);
FluxLine2 = FluxIndex2*(Rmax(2)-Rmin(2))/(N(2));
jj = 1;
FluxQ2 = zeros(n_steps+1, 1);
CumFluxQ2 =  zeros(n_steps+1, 1);


%Total energy 
E_Total = 5.0 + 0.1*iii;                              %In eV
E_Total = Conv.FromeVtoHartree(E_Total);    %In Hartree

disp(['Etotal: ', num2str(E_Total)]);

E_Kin   = E_Total - E_Vib0;            %In Hartree
disp(['Ekin: ', num2str(E_Kin)]);
if E_Kin < 0 
    newData = [Conv.FromHartreetoeV(E_Total)', 0', (iii)'];
    existingData = [existingData; newData];
    continue;
end
k0      = sqrt(2*E_Kin*mu2(1));


psi_0 = exp(1i*k0.* (R1 -x0(1))) .* ...
        exp(-1/2 *((R1  -x0(1))/sigma0(1)).^2) ...
        .*Wavefunction;
% .* exp(-1/2 *((R2-x0(2))/sigma0(2)).^2);

% Normalize
dens0 = psi_0.*conj(psi_0);
N0 = sqrt(trapz(r1, ...
            trapz(r2,dens0,2), ...
          1)  );
psi_0 = psi_0/N0;

% Wavepacket at other times
psi_t = zeros([size(psi_0), n_outSteps+1]);
psi_t(:,:,1) = psi_0;

%% Initialisation of CAPS


% Define the rectangle coordinates
x1ReactantCAP = [7/10, 7/10 1, 1] * Rmax(1);
y1ReactantCAP = [0, 4, 4, 0];
x2ProductCAP = [0, 0, 4/4, 4/4] * Rmax(1);
y2ProductCAP = [6/10, 1, 1, 6/10] * Rmax(2);

%CAPs
L_React = x1ReactantCAP(3) - x1ReactantCAP(2);
L_Prod = y2ProductCAP(2) - y2ProductCAP(1);
React_n   = 4;
Prod_n    = 4;
W_React_L   = (L_React-x1ReactantCAP(1))^(React_n);
W_Prod_L    = (L_Prod-y2ProductCAP(1))^(Prod_n);
alpha_React = factorial(React_n)/(4*1i)*(1i/2)^(React_n);
alpha_Prod = factorial(Prod_n)/(4*1i) * (1i/2)^(Prod_n);
b_React = k0*L_React;
b_Prod = k0*L_Prod;
a_React = lambertw(b_React^(2*React_n+2)/((React_n+1)^(2)*2*alpha_React^(2)))*(React_n+1)/b_React;
a_Prod  = lambertw(b_Prod^(2*Prod_n+2)/((Prod_n+1)^(2)*2*alpha_Prod^(2)))*(Prod_n+1)/b_Prod;
eta_Reactant     = E_Kin*a_React/W_React_L;
eta_Product      = E_Kin*a_Prod/W_Prod_L;

W_React   = (R1-x1ReactantCAP(1)).^(React_n).*heaviside(R1-x1ReactantCAP(1)).*heaviside(y1ReactantCAP(2)-R2);
W_Prod    = (R2-y2ProductCAP(1)).^(Prod_n).*heaviside(x2ProductCAP(3)-R1).*heaviside(R2-y2ProductCAP(1));
%eta_Reactant     = 0.0001;
%eta_Product      = 1;
W_Reactant = eta_Reactant.*W_React;
W_Product = eta_Product.*W_Prod;


%% Operators

dt = tstep * 1E-15 / SI.ta;   % In a.u.

SHH1 = 0.02;
SHH2 = 0.004;
SHH3 = -0.01;


Ebond1 = @(R1,R2)          DeHF*(exp(-2*alphaHF*(R1 - R2*mu(3)/(mu(2)+mu(3)) - ReqHF)) - 2*exp(-alphaHF*(R1 - R2*mu(3)/(mu(2)+mu(3)) - ReqHF)));
Eanti1 = @(R1,R2)      1/2*DeHF*(exp(-2*alphaHF*(R1 - R2*mu(3)/(mu(2)+mu(3)) - ReqHF)) + 2*exp(-alphaHF*(R1 - R2*mu(3)/(mu(2)+mu(3)) - ReqHF)));
Ebond2 = @(R2)          DeHCl*(exp(-2*alphaHCl*(R2 - ReqHCl)) - 2*exp(-alphaHCl*(R2 - ReqHCl)));
Eanti2 = @(R2)      1/2*DeHCl*(exp(-2*alphaHCl*(R2 - ReqHCl)) + 2*exp(-alphaHCl*(R2 - ReqHCl)));
Ebond3 = @(R1, R2)      DeFCl*(exp(-2*alphaFCl*(R1 - R2*mu(3)/(mu(2)+mu(3)) + R2 - ReqFCl)) - 2*exp(-alphaFCl*(R1 - R2*mu(3)/(mu(2)+mu(3)) + R2 - ReqFCl)));
Eanti3 = @(R1, R2)  1/2*DeFCl*(exp(-2*alphaFCl*(R1 - R2*mu(3)/(mu(2)+mu(3)) + R2 - ReqFCl)) + 2*exp(-alphaFCl*(R1 - R2*mu(3)/(mu(2)+mu(3)) + R2 - ReqFCl)));

JAB = @(R1,R2) 1/2*Ebond1(R1,R2) + 1/2*(1 - SHH1)*Eanti1(R1,R2)/(1 + SHH1);
KAB = @(R1,R2) 1/2*Ebond1(R1,R2) - 1/2*(1 - SHH1)*Eanti1(R1,R2)/(1 + SHH1);
JBC = @(R2) 1/2*Ebond2(R2) + 1/2*(1 - SHH2)*Eanti2(R2)/(1 + SHH2);
KBC = @(R2) 1/2*Ebond2(R2) - 1/2*(1 - SHH2)*Eanti2(R2)/(1 + SHH2);
JAC = @(R1, R2) 1/2*Ebond3(R1, R2) + 1/2*(1 - SHH3)*Eanti3(R1, R2)/(1 + SHH3);
KAC = @(R1, R2) 1/2*Ebond3(R1, R2) - 1/2*(1 - SHH3)*Eanti3(R1, R2)/(1 + SHH3);

W = JAB(10000,ReqHCl) + JBC(ReqHCl) + JAC(10000, ReqHCl) - sqrt(1/2*((KAB(10000,ReqHCl) - KBC(ReqHCl))^2 + (KBC(ReqHCl) - KAC(10000, ReqHCl))^2 + (KAC(10000, ReqHCl) - KAB(10000,ReqHCl))^2));
V_First = @(R1, R2) JAB(R1,R2) + JBC(R2) + JAC(R1, R2) - sqrt(1/2*(KAB(R1,R2) - KBC(R2))^2 + 1/2*(KBC(R2) - KAC(R1, R2))^2 + 1/2*(KAC(R1, R2) - KAB(R1,R2))^2) -W;


for i = 1:N(1)
    for j = 1:N(2)
        VContours(i, j) = V_First(R1(i, j), R2(i, j));
    end
end


% *** kinetic energy
T = 0.5 * (P1.^2/mu2(1) + P2.^2/mu2(2));

% Propagator
expV_2  = exp(- 0.5i * VContours * dt/2);
expV    = expV_2.^2                     ;
expT    = exp(- 1i * T * dt/2)          ;
exp_WReactant = exp(-W_Reactant * dt)   ;
exp_WProduct = exp(-W_Product * dt)   ;

%% Propagation
for s = 2:n_outSteps+1
    
    % First potential part
    psi_t(:,:,s) = expV_2 .* psi_t(:,:,s-1);

    for i=1:n_silentSteps-1

        psi_t(:,:,s) = fftn(psi_t(:,:,s));

        %Calculate flux
        flux = psi_t(:,:,s);
        flux = P2 .* flux;
        flux = ifftn(flux)/mu2(2);
        flux = conj(ifftn(psi_t(:,:,s))) .* flux;
        flux = real(flux);
        FluxQ2(jj) = trapz(r1,(flux(:,FluxIndex2)));
        jj = jj + 1;


        %Let's get back to business, I don't got no
        %time to play around, what is this
        psi_t(:,:,s) = expT .* psi_t(:,:,s);
        psi_t(:,:,s) = ifftn(psi_t(:,:,s));
        psi_t(:,:,s) = expV_2 .* psi_t(:,:,s);
        
        

        psi_t(:,:,s) = exp_WProduct .* psi_t(:,:,s);
        psi_t(:,:,s) = exp_WReactant .* psi_t(:,:,s);

        
        psi_t(:,:,s) = expV_2 .* psi_t(:,:,s);
        psi_t(:,:,s) = fftn(psi_t(:,:,s));
        psi_t(:,:,s) = expT .* psi_t(:,:,s);
        psi_t(:,:,s) = ifftn(psi_t(:,:,s));
        psi_t(:,:,s) = expV .* psi_t(:,:,s);
        
    end
    
    % Kinetic part
    psi_t(:,:,s) = fftn(psi_t(:,:,s));
    
    %Calculate flux
    flux = psi_t(:,:,s);
    flux = P2 .* flux;
    flux = ifftn(flux)/mu2(2);
    flux = conj(ifftn(psi_t(:,:,s))) .* flux;
    flux = real(flux);
    FluxQ2(jj) = trapz(r1,(flux(:,FluxIndex2)));
    jj = jj + 1;
    
    
        
    psi_t(:,:,s) = expT .* psi_t(:,:,s);
    psi_t(:,:,s) = ifftn(psi_t(:,:,s));

    % Second potential part
    psi_t(:,:,s) = expV_2 .* psi_t(:,:,s);
    
    
    % Cap part
    psi_t(:,:,s) = exp_WProduct .* psi_t(:,:,s);
    psi_t(:,:,s) = exp_WReactant .* psi_t(:,:,s);
    

    % Third potential part
    psi_t(:,:,s) = expV_2 .* psi_t(:,:,s);
    
    % Second Kinetic part
    psi_t(:,:,s) = fftn(psi_t(:,:,s));
    psi_t(:,:,s) = expT .* psi_t(:,:,s);
    psi_t(:,:,s) = ifftn(psi_t(:,:,s));

    % Fourth potential part
    psi_t(:,:,s) = expV_2 .* psi_t(:,:,s);

    disp( ['Step ', num2str((s-1)*n_silentSteps), ...
           ' out of ', num2str(n_steps), ...
           ' completed' ] );
    
end

%% Show wavepacket



% Axes
t = 0:n_silentSteps*tstep:tfinal;
t_flux = 0:tstep:tfinal;


R1_aa = R1;                 % In angstrom
R2_aa = R2;                 % In angstrom

Contours = Conv.FromeVtoHartree([0.1,0.2,0.4,0.5,0.55,0.57,0.58,0.6,0.7,1,1.25,1.5,2,2.5]);

% Save snapshots and create video
fig_snap = figure;
set(fig_snap, 'NumberTitle', 'off');
set(fig_snap, 'Name', 'Wavepacket Snapshots');
ax_snap = gca;

% Preallocate a video writer object
video_filename = 'wavepacket_video.mp4';
video_writer = VideoWriter(video_filename, 'MPEG-4');
video_writer.FrameRate = 10; % Set the frame rate as needed
open(video_writer);

for t_idx = 1:n_outSteps
    % Plot the wavefunction
    contourf(ax_snap, R1_aa(:,:,1), R2_aa(:,:,1), abs(psi_t(:,:,t_idx)).^2, ...
        [0.001, 0.01, 0.05, 0.1, 0.125, 0.15, 0.2, 0.3, 0.4]);
    
    xlim(ax_snap, [Rmin(1), Rmax(1)]);
    ylim(ax_snap, [Rmin(2), Rmax(2)]);
    
    c = colorbar(ax_snap);
    colormap(ax_snap, 'jet');
    c.Label.String = '$|\psi|^2$';
    c.Label.Interpreter = 'latex';
    ylabel(ax_snap, '$R_{AB}$ [a$_{0}$]');
    xlabel(ax_snap, '$R_{BC}$ [a$_{0}$]');
    title(ax_snap, ['$t = ', num2str(t(t_idx)), '$ fs']);
    hold(ax_snap, 'on');
    
    contour(ax_snap, R1, R2, VContours, Contours);
   
    % Capture the frame for the video
    frame = getframe(fig_snap);
    writeVideo(video_writer, frame);
    hold(ax_snap, 'off');
end

% Close the video writer
close(video_writer);
disp(['Video saved as "', video_filename, '"']);

% Close the snapshot figure
close(fig_snap);


%% Calculate expectation values

FluxQ2(:) = Conv.FromFstohbarEhart(FluxQ2(:));
CumFluxQ2(:) =cumtrapz(t_flux, FluxQ2(:));

newData = [Conv.FromHartreetoeV(E_Total)', max(CumFluxQ2(:))', (iii)'];
existingData = [existingData; newData];


end

% Save the updated data to the text file outside the loop
writematrix(existingData, 'Vib1FHCl.txt', 'Delimiter', '\t');


%% Reset to factory settings

set(groot, 'defaultTextInterpreter', 'remove');
set(groot, 'defaultAxesTickLabelInterpreter','remove');
set(groot, 'defaultColorbarTickLabelInterpreter','remove');
set(groot, 'defaultLegendInterpreter','remove');
set(groot,'defaultAxesFontName','remove');
set(groot,'defaultAxesFontSize', 'remove');