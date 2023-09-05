function plotcase1(z,tvec,RMS12,RMS13,RMS23,orbit_type,ToF,SNC) %3sats

global save_plot

if z == 1
    figure('name','Standard Deviations of DGPS data','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
elseif z == 2
    figure('name','Standard Deviations of ISLS data','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
elseif z == 3
    figure('name','Standard Deviations of DGPS+ISLs data','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
end

hold on; grid on;
subplot(3,1,1)
hold on;
%semilogy(tvec,RMS12,'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,RMS12,'r',tvec,-RMS12,'r','LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,1/3*RMS12,'b--',tvec,-1/3*RMS12,'b--','LineWidth',1); % RMS = 1* sqrt(sigmas)

xlabel('time (s)');
ylabel('standard deviation (m)');
title('RMS SAT1-SAT2');
lgd = legend('3\sigma','-3\sigma','\sigma','-\sigma');
lgd.Interpreter = 'tex';

subplot(3,1,2)
hold on;
%semilogy(tvec,RMS13,'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,RMS13,'r',tvec,-RMS13,'r','LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,1/3*RMS13,'b--',tvec,-1/3*RMS13,'b--', 'LineWidth',1); % RMS = 1* sqrt(sigmas)
xlabel('time (s)');
ylabel('standard deviation (m)');
title('RMS SAT1-SAT3');
lgd = legend('3\sigma','-3\sigma','\sigma','-\sigma');
lgd.Interpreter = 'tex';

subplot(3,1,3)
hold on;
%semilogy(tvec,RMS23,'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,RMS23,'r',tvec,-RMS23,'r','LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,1/3*RMS23,'b--',tvec,-1/3*RMS23,'b--', 'LineWidth',2); % RMS = 1* sqrt(sigmas)
xlabel('time (s)');
ylabel('standard deviation (m)');
title('RMS SAT2-SAT3');
lgd = legend('3\sigma','-3\sigma','\sigma','-\sigma');
lgd.Interpreter = 'tex';

if save_plot == 1
    if z == 1
        if strcmpi(orbit_type,'circular') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/1orbitalperiod/','DGPS3sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/24hours/','DGPS3sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/10days/','DGPS3sats.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','DGPS3sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','DGPS3sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','DGPS3sats.png'));
                end
            end

        elseif strcmpi(orbit_type,'oblate') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/1orbitalperiod/','DGPS3sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/24hours/','DGPS3sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/10days/','DGPS3sats.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','DGPS3sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','DGPS3sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','DGPS3sats.png'));
                end
            end
        end

    elseif z == 2

        if strcmpi(orbit_type,'circular') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/1orbitalperiod/','ISLS3sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/24hours/','ISLS3sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/10days/','ISLS3sats.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','ISLS3sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','ISLS3sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','ISLS3sats.png'));
                end
            end

        elseif strcmpi(orbit_type,'oblate') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/1orbitalperiod/','ISLS3sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/24hours/','ISLS3sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/10days/','ISLS3sats.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','ISLS3sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','ISLS3sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','ISLS3sats.png'));
                end
            end
        end

    elseif z == 3

        if strcmpi(orbit_type,'circular') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/1orbitalperiod/','DGPSISLS3sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/24hours/','DGPSISLS3sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/10days/','DGPSISLS3sats.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','DGPSISLS3sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','DGPSISLS3sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','DGPSISLS3sats.png'));
                end
            end

        elseif strcmpi(orbit_type,'oblate') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/1orbitalperiod/','DGPSISLS3sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/24hours/','DGPSISLS3sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/10days/','DGPSISLS3sats.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','DGPSISLS3sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','DGPSISLS3sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','DGPSISLS3sats.png'));
                end
            end
        end
    end
end
