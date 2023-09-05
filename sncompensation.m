function [Pkbar] = sncompensation(SNC,tt,STMkk_1,Pk_1bar,configuration)
                %% SNC - STATE NOISE COMPENSATION
                %process noise covariance, matrix Q.
                %xdot = Ax + Bu + Te
                %T = process noise mapping matrix
                %e = epsilon = continuos time process noise

                %% table of records
                %sigma = 1e-10 --> final_lambda = 0.1148m,
                %final_frequency = 2.613 GHz
                %
                %sigma = 1e-12 --> final_lambda = 0.1144m,
                %final_frequency = 2.621GHz
                %
                %sigma<1e-12 not worth. 

                if SNC == 1
                    sigma = 1e-10; %[]
                    Q = sigma^2*eye(3); %3x3

                    if (tt ==1)
                        T = (tt)* [(tt)/2 * eye(3); eye(3)]; %6x3
                    else
                        T = (tt - (tt-1))* [(tt - (tt-1))/2 * eye(3); eye(3)]; %6x3
                    end

                    g = T*Q*T'; %6x6
                    if configuration == 1
                        gamma = blkdiag(g,g,g); %18x18, 3 sats
                    elseif configuration == 2
                        gamma = blkdiag(g,g,g,g); %24x24, 4 sats
                    elseif configuration == 3
                        gamma = blkdiag(g,g,g,g,g); %30x30, 5 sats
                    end
                    Pkbar   = STMkk_1 * Pk_1bar * STMkk_1' + gamma; %18x18, Propagate a-priori state deviation vector and Covariance Matrix

                else
                    Pkbar   = STMkk_1 * Pk_1bar * STMkk_1'; %18x18, Propagate a-priori state deviation vector and Covariance Matrix
                end

        