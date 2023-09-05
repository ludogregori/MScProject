function plotcase3(z,tvec,RMS10,RMS20,RMS30,RMS40,orbit_type,ToF,SNC) %5sats: 4sats + 1central

global save_plot

if     z == 1
    figure('name','Standard Deviations of DGPS data','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
elseif z == 2
    figure('name','Standard Deviations of ISLS data','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
elseif z == 3
    figure('name','Standard Deviations of DGPS+ISLs data','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
end

hold on; grid on;

subplot(1,4,1)
hold on;
%semilogy(tvec,RMS10,'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,RMS10,'r',tvec,-RMS10,'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,1/3*RMS10,'b--',tvec,-1/3*RMS10,'b--','LineWidth',1); % RMS = 1* sqrt(sigmas)
title('SAT1 - SAT2');
xlabel('time (s)');
ylabel('standard deviation (m)');
%lgd = legend('3\sigma of DGPS','3\sigma of DGPS+ISLs');
lgd = legend('3\sigma of DGPS','-3\sigma of DGPS+ISLs','1\sigma of DGPS','-1\sigma of DGPS+ISLs');
lgd.Interpreter = 'tex';

subplot(1,4,2)
hold on;
%semilogy(tvec,RMS20,'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,RMS20,'r',tvec,-RMS20,'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,1/3*RMS20,'b--',tvec,-1/3*RMS20,'b--','LineWidth',1); % RMS = 1* sqrt(sigmas)
title('SAT2 - SAT3');
xlabel('time (s)');
ylabel('standard deviation (m)');
%lgd = legend('3\sigma of DGPS','3\sigma of DGPS+ISLs');
lgd = legend('3\sigma of DGPS','-3\sigma of DGPS+ISLs','1\sigma of DGPS','-1\sigma of DGPS+ISLs');
lgd.Interpreter = 'tex';

subplot(1,4,3)
hold on;
%semilogy(tvec,RMS30,'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,RMS30,'r',tvec,-RMS30,'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,1/3*RMS30,'b--',tvec,-1/3*RMS30,'b--','LineWidth',1); % RMS = 1* sqrt(sigmas)
title('SAT3 - SAT4');
xlabel('time (s)');
ylabel('standard deviation (m)');
%lgd = legend('3\sigma of DGPS','3\sigma of DGPS+ISLs');
lgd = legend('3\sigma of DGPS','-3\sigma of DGPS+ISLs','1\sigma of DGPS','-1\sigma of DGPS+ISLs');
lgd.Interpreter = 'tex';

subplot(1,4,4)
hold on;
%semilogy(tvec,RMS40,'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,RMS40,'r',tvec,-RMS40,'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
plot(tvec,1/3*RMS40,'b--',tvec,-1/3*RMS40,'b--','LineWidth',1); % RMS = 1* sqrt(sigmas)
title('SAT4 - SAT1');
xlabel('time (s)');
ylabel('standard deviation (m)');
%lgd = legend('3\sigma of DGPS','3\sigma of DGPS+ISLs');
lgd = legend('3\sigma of DGPS','-3\sigma of DGPS+ISLs','1\sigma of DGPS','-1\sigma of DGPS+ISLs');
lgd.Interpreter = 'tex';

hold off;

if save_plot == 1
    if z == 1
        if strcmpi(orbit_type,'circular') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/1orbitalperiod/','DGPS5sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/24hours/','DGPS5sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/10days/','DGPS5sats.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/1orbitalperiod/','DGPS5sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/24hours/','DGPS5sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/10days/','DGPS5sats.png'));
                end
            end

        elseif strcmpi(orbit_type,'oblate') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config3/1orbitalperiod/','DGPS5sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config3/24hours/','DGPS5sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config3/10days/','DGPS5sats.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/1orbitalperiod/','DGPS5sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/24hours/','DGPS5sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/10days/','DGPS5sats.png'));
                end
            end
        end

    elseif z == 2

        if strcmpi(orbit_type,'circular') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/1orbitalperiod/','ISLS5sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/24hours/','ISLS5sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/10days/','ISLS5sats.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/1orbitalperiod/','ISLS5sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/24hours/','ISLS5sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/10days/','ISLS5sats.png'));
                end
            end

        elseif strcmpi(orbit_type,'oblate') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config3/1orbitalperiod/','ISLS5sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config3/24hours/','ISLS5sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config3/10days/','ISLS5sats.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/1orbitalperiod/','ISLS5sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/24hours/','ISLS5sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/10days/','ISLS5sats.png'));
                end
            end
        end

    elseif z == 3

        if strcmpi(orbit_type,'circular') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/1orbitalperiod/','DGPSISLS5sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/24hours/','DGPSISLS5sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/10days/','DGPSISLS5sats.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/1orbitalperiod/','DGPSISLS5sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/24hours/','DGPSISLS5sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/10days/','DGPSISLS5sats.png'));
                end
            end

        elseif strcmpi(orbit_type,'oblate') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config3/1orbitalperiod/','DGPSISLS5sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config3/24hours/','DGPSISLS5sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config3/10days/','DGPSISLS5sats.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/1orbitalperiod/','DGPSISLS5sats.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/24hours/','DGPSISLS5sats.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/snc ON/10days/','DGPSISLS5sats.png'));
                end
            end
        end
    end
end