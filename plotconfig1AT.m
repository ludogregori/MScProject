function plotconfig1AT(tvec,RMS,orbit_type,ToF,SNC)

global save_plot

figure('name','Standard Deviations of DGPS vs ISLS vs DGPS+ISLs data (log scale)','units','normalized','outerposition',[0.2 0.2 0.8 0.8]); %figure 7
hold on; grid on;

subplot(3,1,1)
semilogy(tvec,abs(RMS{1}(:,1)),tvec,abs(RMS{1}(:,2)),tvec,abs(RMS{1}(:,3)),'r', 'LineWidth',2); % RMS = 3*sqrt(sigmas)
title('DGPS');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of DGPS');
lgd.Interpreter = 'tex';

subplot(3,1,2)
semilogy(tvec,abs(RMS{2}(:,1)),tvec,abs(RMS{2}(:,2)),tvec,abs(RMS{2}(:,3)),'b', 'LineWidth',2); % RMS = 3*sqrt(sigmas)
title('ISLS');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of ISLS');
lgd.Interpreter = 'tex';

subplot(3,1,3)
semilogy(tvec,abs(RMS{3}(:,1)),tvec,abs(RMS{3}(:,2)),tvec,abs(RMS{3}(:,3)),'k', 'LineWidth',2); % RMS = 3*sqrt(sigmas)
title('DGPS + ISLS');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of DGPS+ISLS');
lgd.Interpreter = 'tex';

if save_plot == 1
    if strcmpi(orbit_type,'circular') == 1
        if SNC == 0
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/1orbitalperiod/','allsystems3sat.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/24hours/','allsystems3sat.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/10days/','allsystems3sat.png'));
            end
        elseif SNC == 1
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','allsystems3sat.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','allsystems3sat.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','allsystems3sat.png'));
            end
        end

    elseif strcmpi(orbit_type,'oblate') == 1
        if SNC == 0
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/1orbitalperiod/','allsystems3sat.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/24hours/','allsystems3sat.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/10days/','allsystems3sat.png'));
            end
        elseif SNC == 1
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','allsystems3sat.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','allsystems3sat.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','allsystems3sat.png'));
            end
        end
    end
end
