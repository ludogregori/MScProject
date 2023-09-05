function [STMkk_1] = STMatrix(tvec,nc,tt,configuration)

c = cos(nc*tvec(tt));
s = sin(nc*tvec(tt));
STM = [4-3*c             0    0      s/nc           (2/nc)*(1-c)      0;
    6*(s-nc*tvec(tt)) 1    0 -2*(1-c)/nc   (4*s-3*nc*tvec(tt))/nc  0;
    0                 0    c       0               0             s/nc;
    3*nc*s            0    0       c              2*s              0;
    -6*nc*(1-c)        0    0     -2*s            4*c-3             0;
    0                 0  -nc*s     0               0               c]; %STM matrix

if configuration == 1
    newSTM(:,:,tt) = blkdiag(STM,STM,STM); %18x18, 3 sats
elseif configuration == 2
    newSTM(:,:,tt) = blkdiag(STM,STM,STM,STM); %24x24, 4 sats
elseif configuration == 3
    newSTM(:,:,tt) = blkdiag(STM,STM,STM,STM,STM); %30x30, 5 sats
end
if configuration == 1
    newSTM(:,:,tt) = blkdiag(STM,STM,STM); %18x18

elseif configuration == 2
    newSTM(:,:,tt) = blkdiag(STM,STM,STM,STM); %24x24

elseif configuration == 3
    newSTM(:,:,tt) = blkdiag(STM,STM,STM,STM,STM); %30x30
end

cellSTM(tt) = {newSTM(:,:,tt)}; %storing the STM in a cell array

if(tt == 1)
    STMkk_1 = newSTM(:,:,tt);
else
    STMkk_1 = newSTM(:,:,tt) / newSTM(:,:,tt-1);  % STM(tk,tk_1) = STM(tk, 0) / STM(tk_1, 0), conventional KF
end