function plotXALLconfig1(z,tvec,RMS,orbit_type,ToF,SNC)
%config1 , 3 sats
global save_plot

if z == 1
    figure('name','Standard Deviations of DGPS data','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
elseif z == 2
    figure('name','Standard Deviations of ISLS data','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
elseif z == 3
    figure('name','Standard Deviations of DGPS+ISLs data','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
end

%for i = 1 : length(tvec)

hold on;
subplot(1,3,1) %DPGS,ISLS,DGPS + ISLS, 1 sat
hold on;grid on;grid minor;
loglog(tvec,RMS{1}(:,1)','r',tvec,RMS{2}(:,1)','b',tvec,RMS{3}(:,1)','g','LineWidth',2); % RMS = 3* sqrt(sigmas)
%loglog(tvec,1/3*RMS{1}(:,1)','r',tvec,1/3*RMS{2}(:,1)','b',tvec,1/3*RMS{3}(:,1)','g','LineWidth',2); % RMS = 1* sqrt(sigmas)
xlabel('time (s)');
ylabel('standard deviation (m)');
xticks([5700 86400 864000]);
xline([5700 86400 864000],'-k',{'1 orbital day','24 hours','10 days'},'Linewidth',2);
xl.LabelHorizontalAlignment = 'center';
xticklabels({'1 orbital day','24 hours','10 days'})
title('RMS SAT1-SAT2');
lgd = legend('3\sigma DGPS','3\sigma ISLS','3\sigma DGPS + ISLS','Location','best');
lgd.Interpreter = 'tex';
set(gca, 'XScale', 'log', 'YScale', 'log');

subplot(1,3,2) %DPGS,ISLS,DGPS + ISLS 2 sat
hold on;grid on;grid minor;
%semilogy(tvec,RMS13,'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
loglog(tvec,RMS{1}(:,2)','r',tvec,RMS{2}(:,2)','b',tvec,RMS{3}(:,2)','g','LineWidth',2); % RMS = 3* sqrt(sigmas)
%loglog(tvec,1/3*RMS{1}(:,2)','r',tvec,1/3*RMS{2}(:,2)','b',tvec,1/3*RMS{3}(:,2)','g','LineWidth',2); % RMS = 1* sqrt(sigmas)
xlabel('time (s)');
ylabel('standard deviation (m)');
xticks([5700 86400 864000]);
xline([5700 86400 864000],'-k',{'1 orbital day','24 hours','10 days'},'Linewidth',2);
xticklabels({'1 orbital day','24 hours','10 days'})
title('RMS SAT1-SAT3');
lgd = legend('3\sigma DGPS','3\sigma ISLS','3\sigma DGPS + ISLS','Location','best');
lgd.Interpreter = 'tex';
set(gca, 'XScale', 'log', 'YScale', 'log');

subplot(1,3,3) %DPGS,ISLS,DGPS + ISLS 3 sat
hold on;grid on;grid minor;
%semilogy(tvec,RMS23,'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
loglog(tvec,RMS{1}(:,3)','r',tvec,RMS{2}(:,3)','b',tvec,RMS{3}(:,3)','g','LineWidth',2); % RMS = 3* sqrt(sigmas)
%loglog(tvec,1/3*RMS{1}(:,3)','r',tvec,1/3*RMS{2}(:,3)','b',tvec,1/3*RMS{3}(:,3)','g','LineWidth',2); % RMS = 1* sqrt(sigmas)
xlabel('time (s)');
ylabel('standard deviation (m)');
xticks([5700 86400 864000]);
xline([5700 86400 864000],'-k',{'1 orbital day','24 hours','10 days'},'Linewidth',2);
xticklabels({'1 orbital day','24 hours','10 days'})
title('RMS SAT2-SAT3');
lgd = legend('3\sigma DGPS','3\sigma ISLS','3\sigma DGPS + ISLS','Location','best');
lgd.Interpreter = 'tex';
set(gca, 'XScale', 'log', 'YScale', 'log');
%end

if save_plot == 1
    if z == 1
        if strcmpi(orbit_type,'circular') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/1orbitalperiod/','satcomparison.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/24hours/','satcomparison.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/10days/','satcomparison.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','satcomparison.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','satcomparison.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','satcomparison.png'));
                end
            end

        elseif strcmpi(orbit_type,'oblate') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/1orbitalperiod/','satcomparison.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/24hours/','satcomparison.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/10days/','satcomparison.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','satcomparison.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','satcomparison.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','satcomparison.png'));
                end
            end
        end

    elseif z == 2

        if strcmpi(orbit_type,'circular') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/1orbitalperiod/','satcomparison.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/24hours/','satcomparison.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/10days/','satcomparison.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','satcomparison.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','satcomparison.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','satcomparison.png'));
                end
            end

        elseif strcmpi(orbit_type,'oblate') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/1orbitalperiod/','satcomparison.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/24hours/','satcomparison.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/10days/','satcomparison.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','satcomparison.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','satcomparison.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','satcomparison.png'));
                end
            end
        end

    elseif z == 3

        if strcmpi(orbit_type,'circular') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/1orbitalperiod/','satcomparison.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/24hours/','satcomparison.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/10days/','satcomparison.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','satcomparison.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','satcomparison.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','satcomparison.png'));
                end
            end

        elseif strcmpi(orbit_type,'oblate') == 1
            if SNC == 0
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/1orbitalperiod/','satcomparison.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/24hours/','satcomparison.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/10days/','satcomparison.png'));
                end
            elseif SNC == 1
                if ToF == 1
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','satcomparison.png'));
                elseif ToF == 24
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','satcomparison.png'));
                elseif ToF == 10
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','satcomparison.png'));
                end
            end
        end
    end
end
