syms R1 R2

f = @(R1,R2) V(R1, R2) - W;
Saddle_Bool = 1;
Z = V(R1, R2) - W;
grad = gradient(Z, [R1, R2]);
for n = 1:10
    RCrit = vpasolve(grad == 0, [R1 R2], RangeGuess);%, 'Random', true);
    
    SU = subs(det(hessian(Z, [R1, R2])), ...
        [R1, R2], [RCrit.R1, RCrit.R2]);
    
    if SU < 0
        if RCrit.R1 < RCrit.R2
            fprintf("Late barrier \n")
        end
        if RCrit.R1 > RCrit.R2
            fprintf("Early barrier \n")
        end
        fprintf("The saddle point has the coordinates " + ...
        "provided by: \n%f [Å] and %f [Å] \n", RCrit.R1, RCrit.R2);
        fprintf("The barrier height is: \n%f [eV] \n", f(RCrit.R1, RCrit.R2));
        Saddle_Bool = 0;
        break;
    end
    if SU >= 0
        Saddle_Bool = 2;
        break;
    end
end

if Saddle_Bool == 1
    disp("No critical point");
    return;
end 

if Saddle_Bool ~= 2 
    return
end
MaxMin = subs(grad(1), {R1,R2},{RCrit.R1, RCrit.R2});

if  MaxMin >= 0
    fprintf("You found the minimum coordinates: \n%f and %f\n", R1Crit, R2Crit);
    return
end
if  MaxMin <= 0
    fprintf("You found the maximum coordinates: \n%f %f\n", R1Crit, R2Crit);
    return
end
if Answer == 0 
    fprintf("Inconclusive results");
    return
end
