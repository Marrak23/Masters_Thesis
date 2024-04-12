Parameters
SHH1 = 0.1;
SHH2 = 0.1;
SHH3 = 0.1;



Ebond1 = @(R1) DeHH*(exp(-2*alphaHH*(R1 - ReqHH)) - 2*exp(-alphaHH*(R1 - ReqHH)));
Eanti1 = @(R1) 1/2*DeHH*(exp(-2*alphaHH*(R1 - ReqHH)) + 2*exp(-alphaHH*(R1 - ReqHH)));
Ebond2 = @(R2) DeHH*(exp(-2*alphaHH*(R2 - ReqHH)) - 2*exp(-alphaHH*(R2 - ReqHH)));
Eanti2 = @(R2) 1/2*DeHH*(exp(-2*alphaHH*(R2 - ReqHH)) + 2*exp(-alphaHH*(R2 - ReqHH)));
Ebond3 = @(R1, R2) DeHH*(exp(-2*alphaHH*(R1 + R2 - ReqHH)) - 2*exp(-alphaHH*(R1 + R2 - ReqHH)));
Eanti3 = @(R1, R2) 1/2*DeHH*(exp(-2*alphaHH*(R1 + R2 - ReqHH)) + 2*exp(-alphaHH*(R1 + R2 - ReqHH)));


JAB = @(R1) 1/2*Ebond1(R1) + 1/2*(1 - SHH1)*Eanti1(R1)/(1 + SHH1);
KAB = @(R1) 1/2*Ebond1(R1) - 1/2*(1 - SHH1)*Eanti1(R1)/(1 + SHH1);
JBC = @(R2) 1/2*Ebond2(R2) + 1/2*(1 - SHH2)*Eanti2(R2)/(1 + SHH2);
KBC = @(R2) 1/2*Ebond2(R2) - 1/2*(1 - SHH2)*Eanti2(R2)/(1 + SHH2);
JAC = @(R1, R2) 1/2*Ebond3(R1, R2) + 1/2*(1 - SHH3)*Eanti3(R1, R2)/(1 + SHH3);
KAC = @(R1, R2) 1/2*Ebond3(R1, R2) - 1/2*(1 - SHH3)*Eanti3(R1, R2)/(1 + SHH3);

V = @(R1, R2) JAB(R1) + JBC(R2) + JAC(R1, R2) - ...
sqrt(1/2*(KAB(R1) - KBC(R2))^2 + 1/2*(KBC(R2) - KAC(R1, R2))^2 + 1/2*(KAC(R1, R2) - KAB(R1))^2);
W = JAB(10000) + JBC(ReqHH) + JAC(10000, ReqHH) - ... 
sqrt(1/2*((KAB(10000) - KBC(ReqHH))^2 + (KBC(ReqHH) - KAC(10000, ReqHH))^2 + (KAC(10000, ReqHH) - KAB(10000))^2));



x = linspace(0.3, 3.2, 500);
y = linspace(0.3, 3.2, 500);
[X, Y] = meshgrid(x, y);
Z = zeros(size(X));close

for i = 1:numel(x)
    for j = 1:numel(y)
        Z(i, j) = V(X(i, j), Y(i, j)) - W;
    end
end


% Coordinates to add a dot
dotX = 0.95;
dotY = 0.95;

A = [0.57,1.5,1.25,1,0.7,0.6,0.58,0.55,0.5,0.4,0.2,0.1];
contour(X, Y, Z, A);
%title('Sato parameters [S$_{1}$, S$_{2}$, S$_{3}$] = [0.1 ; 0.1 ; 0.1]', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('R$_{\mathrm{HH}}, [\mathrm{\AA}]$','interpreter','latex', 'FontSize', 20);
ylabel('R$_{\mathrm{HH}}, [\mathrm{\AA}]$','interpreter','latex', 'FontSize', 20);
c = colorbar;
ylabel(c, 'V$_{\mathrm{LEPS}}$, [eV]', 'Interpreter', 'latex', 'FontSize', 20);


% Hold on to add a dot
hold on;

% Add a dot at the specified coordinates
scatter(dotX, dotY, 100, 'r', 'filled', 'DisplayName', 'Saddle Point');

%Add reaction text
text(1.4, 2.0, 'H + H$_{2}$ $\rightarrow$ $ $H$_{2}$ + H', 'Interpreter', 'latex', ...
    'FontSize', 26, 'VerticalAlignment', 'bottom');

exportgraphics(gcf, 'HHHRel.png','Resolution',1500);

% Hold off to stop overlaying additional elements
hold off;


