function dX = HCWint(tvec,Xhcw,pars,Nsat)

% integrating via ODE the linearised problem to obtain the trajectory

nc = pars.nc;

% A matrix 
% dx = A*x + B*u

A = zeros(6); %initialiasing the matrix
A(1:3,4:6) = eye(3);
A(4,1) = 3*nc^2;
A(4,5) = 2*nc;
A(5,4) = -2*nc;
A(6,3) = -nc^2;

switch Nsat
%integrating the dynamics
    case 3
        dX1 = A*Xhcw(1:6,1);   %first sat
        dX2 = A*Xhcw(7:12,1);  %second sat
        dX3 = A*Xhcw(13:18,1); %third sat
        
        dX = [dX1; dX2; dX3]; %output vector as a column, 18x1 
    case 4
        dX1 = A*Xhcw(1:6,1);   %first sat
        dX2 = A*Xhcw(7:12,1);  %second sat
        dX3 = A*Xhcw(13:18,1); %third sat
        dX4 = A*Xhcw(19:24,1); %fourth sat
        dX = [dX1; dX2; dX3; dX4]; %output vector as a column, 24x1 
end