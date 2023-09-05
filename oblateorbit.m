function [IC] = oblateorbit(Nsat,Rform,nc,tvec,theta0)
    
for i = 1:Nsat
    IC(1,i) = 0.5 * Rform * cos(nc*tvec(1) + theta0(i)); %x
    IC(2,i) = Rform * sin(nc*tvec(1) + theta0(i)); %y
    IC(3,i) = sqrt(3)/2 * Rform * cos(nc*tvec(1) + theta0(i)); %z
    IC(4,i) = nc*Rform/2 * sin(nc*tvec(1) + theta0(i)); %dx
    IC(5,i) = -nc * Rform * cos(nc*tvec(1) + theta0(i)); %dy
    IC(6,i) = sqrt(3)/2 * nc * Rform * sin(nc*tvec(1) + theta0(i)); %dz
    IC(:,i) = [IC(:,i)];
end