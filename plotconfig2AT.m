function plotconfig2AT(tvec,RMS,orbit_type,ToF,SNC)

global save_plot

figure('name','Standard Deviations of DPGS vs ISLS vs DGPS+ISLs data (log scale)','units','normalized','outerposition',[0.2 0.2 0.8 0.8]); %figure 5 for config 1, figure 6 for config 2/3
hold on; grid on;

subplot(3,1,1)
semilogy(tvec,RMS{1}(:,1),tvec,RMS{1}(:,2),tvec,RMS{1}(:,3),tvec,RMS{1}(:,4),'r', 'LineWidth',2); %RMS = 3* sqrt(sigmas)
title('DGPS');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of DGPS');
lgd.Interpreter = 'tex';

subplot(3,1,2)
semilogy(tvec,RMS{2}(:,1),tvec,RMS{2}(:,2),tvec,RMS{2}(:,3),tvec,RMS{2}(:,4),'b', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
title('ISLS');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of ISLS');
lgd.Interpreter = 'tex';

subplot(3,1,3)
semilogy(tvec,RMS{3}(:,1),tvec,RMS{3}(:,2),tvec,RMS{3}(:,3),tvec,RMS{3}(:,4),'k', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
title('DGPS + ISLS');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of DGPS+ISLS');
lgd.Interpreter = 'tex';

if save_plot == 1
    if strcmpi(orbit_type,'circular') == 1
        if SNC == 0
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/1orbitalperiod/','allsystems4sats.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/24hours/','allsystems4sats.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/10days/','allsystems4sats.png'));
            end
        elseif SNC == 1
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/snc ON/1orbitalperiod/','allsystems4sats.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/snc ON/24hours/','allsystems4sats.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/snc ON/10days/','allsystems4sats.png'));
            end
        end

    elseif strcmpi(orbit_type,'oblate') == 1
        if SNC == 0
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config2/1orbitalperiod/','allsystems4sats.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config2/24hours/','allsystems4sats.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config2/10days/','allsystems4sats.png'));
            end
        elseif SNC == 1
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/snc ON/1orbitalperiod/','allsystems4sats.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/snc ON/24hours/','allsystems4sats.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/snc ON/10days/','allsystems4sats.png'));
            end
        end
    end
end
