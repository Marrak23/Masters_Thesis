%Barrier height is 0.5416 eV
%Barrier height is here 0.204875 eV
%Located at r1 = rOH = 1.23 [Å]
%Located at r2 = rHH = 0.92 [Å]
%Reaction O + H2 ---> OH + H

Parameters
ConstantsSI
SHH1 = 0.15;
SHH2 = 0.4;
SHH3 = -0.455;

mu = [muHO,muHH,muHO]/SI.m_e;


Ebond1 = @(R1,R2)          DeHO*(exp(-2*alphaHO*(R1 - R2*mu(3)/(mu(2)+mu(3)) - ReqHO)) - 2*exp(-alphaHO*(R1 - R2*mu(3)/(mu(2)+mu(3)) - ReqHO)));
Eanti1 = @(R1,R2)      1/2*DeHO*(exp(-2*alphaHO*(R1 - R2*mu(3)/(mu(2)+mu(3)) - ReqHO)) + 2*exp(-alphaHO*(R1 - R2*mu(3)/(mu(2)+mu(3)) - ReqHO)));
Ebond2 = @(R2)          DeHH*(exp(-2*alphaHH*(R2 - ReqHH)) - 2*exp(-alphaHH*(R2 - ReqHH)));
Eanti2 = @(R2)      1/2*DeHH*(exp(-2*alphaHH*(R2 - ReqHH)) + 2*exp(-alphaHH*(R2 - ReqHH)));
Ebond3 = @(R1, R2)      DeHO*(exp(-2*alphaHO*(R1 - R2*mu(3)/(mu(2)+mu(3))  + R2 - ReqHO)) - 2*exp(-alphaHO*(R1 - R2*mu(3)/(mu(2)+mu(3)) + R2 - ReqHO)));
Eanti3 = @(R1, R2)  1/2*DeHO*(exp(-2*alphaHO*(R1 - R2*mu(3)/(mu(2)+mu(3)) + R2 - ReqHO)) + 2*exp(-alphaHO*(R1 - R2*mu(3)/(mu(2)+mu(3)) + R2 - ReqHO)));

JAB = @(R1,R2) 1/2*Ebond1(R1,R2) + 1/2*(1 - SHH1)*Eanti1(R1,R2)/(1 + SHH1);
KAB = @(R1,R2) 1/2*Ebond1(R1,R2) - 1/2*(1 - SHH1)*Eanti1(R1,R2)/(1 + SHH1);
JBC = @(R2) 1/2*Ebond2(R2) + 1/2*(1 - SHH2)*Eanti2(R2)/(1 + SHH2);
KBC = @(R2) 1/2*Ebond2(R2) - 1/2*(1 - SHH2)*Eanti2(R2)/(1 + SHH2);
JAC = @(R1, R2) 1/2*Ebond3(R1, R2) + 1/2*(1 - SHH3)*Eanti3(R1, R2)/(1 + SHH3);
KAC = @(R1, R2) 1/2*Ebond3(R1, R2) - 1/2*(1 - SHH3)*Eanti3(R1, R2)/(1 + SHH3);


W = JAB(10000,ReqHH) + JBC(ReqHH) + JAC(10000, ReqHO) - sqrt(1/2*((KAB(10000,ReqHH) - KBC(ReqHH))^2 + (KBC(ReqHH) - KAC(10000, ReqHH))^2 + (KAC(10000, ReqHH) - KAB(10000,ReqHH))^2));
V = @(R1, R2) JAB(R1,R2) + JBC(R2) + JAC(R1, R2) - sqrt(1/2*(KAB(R1,R2) - KBC(R2))^2 + 1/2*(KBC(R2) - KAC(R1, R2))^2 + 1/2*(KAC(R1, R2) - KAB(R1,R2))^2) - W;


x = linspace(1, 5, 500);
y = linspace(0.25, 3, 500);
[X, Y] = meshgrid(x, y);
Z = zeros(size(X));

for i = 1:numel(x)
    for j = 1:numel(y)
        Z(i, j) = V(X(i, j), Y(i, j));
    end
end

% Coordinates to add a dot
dotX = 1.231371 + 0.920176*mu(3)/(mu(2)+mu(3));
dotY = 0.920176;

A = [1, 0.5,0.45,0.44, 0.4, 0.2, 0.1, 0, -0.5] + 0.2170;
contour(X, Y, Z, A);
%title('Sato parameters [S$_{1}$, S$_{2}$, S$_{3}$] = [0.15 ; 0.4 ; -0.455]', ...
%    'Interpreter', 'latex', 'FontSize', 16);
text(2.75, 2.0, 'O + H$_{2}\rightarrow$ OH + H', 'Interpreter', 'latex', ...
    'FontSize', 24, 'VerticalAlignment', 'bottom');
xlabel('R$_{\mathrm{OH}} [\mathrm{\AA}]$','interpreter','latex', 'FontSize', 20);
ylabel('R$_{\mathrm{HH}} [\mathrm{\AA}]$','interpreter','latex', 'FontSize', 20);
c = colorbar;
ylabel(c, 'V$_{\mathrm{LEPS}}$, [eV]', 'Interpreter', 'latex', 'FontSize', 20);


% Hold on to add a dot
hold on;

% Add a dot at the specified coordinates
scatter(dotX, dotY, 100, 'r', 'filled', 'DisplayName', 'Saddle Point');

%Add line through plot
a = (dotY-min(y))/(dotX-min(x));
b = min(y) - a*min(x);
%line([min(x), (dotX-min(x))*(max(y)-(min(y) - (dotY-min(y))/(dotX-min(x))*min(x)))/(dotY-min(y))],[min(y), max(y)], 'Color','r', 'LineStyle', '--');

%exportgraphics(gcf, 'HHOJac.png','Resolution',1500);

% Hold off to stop overlaying additional elements
hold off;

