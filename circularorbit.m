function [IC] = circularorbit(Nsat,Rform,nc,tvec,theta0)
    
for i = 1:Nsat
    IC(1,i) = 0.5 * Rform * cos(-nc*tvec(1) + theta0(i)); %x
    IC(2,i) =  Rform * sin(-nc*tvec(1) + theta0(i)); %y
    IC(3,i) =  Rform * cos(-nc*tvec(1) + theta0(i)); %z
    IC(4,i) = +0.5 * nc * Rform * sin(-nc*tvec(1) + theta0(i)); %dx
    IC(5,i) = -nc * Rform * cos(-nc*tvec(1) + theta0(i)); %dy
    IC(6,i) = +nc * Rform * sin(-nc*tvec(1) + theta0(i)); %dz
    IC(:,i) = [IC(:,i)];
end