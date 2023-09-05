function plotcase2DGPSvsDGPSISLS(tvec,RMS,orbit_type,ToF,SNC)

global save_plot

figure('name','Standard Deviations of DGPS vs DGPS+ISLS data (log scale)','units','normalized','outerposition',[0.2 0.2 0.8 0.8]); %figure 5 for config 2/3
hold on; grid on;

subplot(6,1,1)
semilogy(tvec,RMS{1}(:,1),tvec,RMS{3}(:,1),'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
title('SAT1 - SAT2');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of ISLS','3\sigma of DGPS+ISLs');
lgd.Interpreter = 'tex';

subplot(6,1,2)
semilogy(tvec,RMS{1}(:,2),tvec,RMS{3}(:,2),'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
title('SAT2 - SAT3');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of DGPS','3\sigma of DGPS+ISLs');
lgd.Interpreter = 'tex';

subplot(6,1,3)
semilogy(tvec,RMS{1}(:,3),tvec,RMS{3}(:,3),'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
title('SAT3 - SAT4');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of DGPS','3\sigma of DGPS+ISLs');
lgd.Interpreter = 'tex';

subplot(6,1,4)
semilogy(tvec,RMS{1}(:,4),tvec,RMS{3}(:,4),'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
title('SAT4 - SAT1');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of DGPS','3\sigma of DGPS+ISLs');
lgd.Interpreter = 'tex';

subplot(6,1,5)
semilogy(tvec,RMS{1}(:,5),tvec,RMS{3}(:,5),'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
title('SAT1 - SAT3');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of DGPS','3\sigma of DGPS+ISLs');
lgd.Interpreter = 'tex';

subplot(6,1,6)
semilogy(tvec,RMS{1}(:,6),tvec,RMS{3}(:,6),'r', 'LineWidth',2); % RMS = 3* sqrt(sigmas)
title('SAT4 - SAT2');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of DGPS','3\sigma of DGPS+ISLs');
lgd.Interpreter = 'tex';

hold off;
if save_plot == 1
    if strcmpi(orbit_type,'circular') == 1
        if SNC == 0
            if ToF == 1
                disp('dddddddddddddd')
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/1orbitalperiod/','DGPSvsDGPSISLS4sats.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/24hours/','DGPSvsDGPSISLS4sats.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/10days/','DGPSvsDGPSISLS4sats.png'));
            end
        elseif SNC == 1
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/snc ON/1orbitalperiod/','DGPSvsDGPSISLS4sats.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/snc ON/24hours/','DGPSvsDGPSISLS4sats.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/snc ON/10days/','DGPSvsDGPSISLS4sats.png'));
            end
        end

    elseif strcmpi(orbit_type,'oblate') == 1
        if SNC == 0
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config2/1orbitalperiod/','DGPSvsDGPSISLS4sats.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config2/24hours/','DGPSvsDGPSISLS4sats.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config2/10days/','DGPSvsDGPSISLS4sats.png'));
            end
        elseif SNC == 1
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/snc ON/1orbitalperiod/','DGPSvsDGPSISLS4sats.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/snc ON/24hours/','DGPSvsDGPSISLS4sats.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/snc ON/10days/','DGPSvsDGPSISLS4sats.png'));
            end
        end
    end
end
