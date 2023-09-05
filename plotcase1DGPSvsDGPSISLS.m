function plotcase1DGPSvsDGPSISLS(tvec,RMS,orbit_type,ToF,SNC)

global save_plot

figure('name','Standard Deviations of DGPS vs DGPS+ISLs data','units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
hold on; grid on;

subplot(3,1,1)
semilogy(tvec,RMS{1}(:,1),tvec,RMS{3}(:,1),'r', 'LineWidth',2); % RMS = 3*sqrt(sigmas)
title('SAT1 - SAT2');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of DGPS','3\sigma of DGPS+ISLs');
lgd.Interpreter = 'tex';

subplot(3,1,2)
semilogy(tvec,RMS{1}(:,2),tvec,RMS{3}(:,2),'r', 'LineWidth',2); % RMS = 3*sqrt(sigmas)
title('SAT1 - SAT3');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of DGPS','3\sigma of DGPS+ISLs');
lgd.Interpreter = 'tex';

subplot(3,1,3)
semilogy(tvec,RMS{1}(:,3),tvec,RMS{3}(:,3),'r', 'LineWidth',2); % RMS = 3*sqrt(sigmas)
title('SAT2 - SAT3');
xlabel('time (s)');
ylabel('standard deviation (m)');
lgd = legend('3\sigma of DGPS','3\sigma of DGPS+ISLs');
lgd.Interpreter = 'tex';

if save_plot == 1
    if strcmpi(orbit_type,'circular') == 1
        if SNC == 0
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/1orbitalperiod/','DGPSvDGPSISLS3sats.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/24hours/','DGPSvDGPSISLS3sats.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/10days/','DGPSvDGPSISLS3sats.png'));
            end
        elseif SNC == 1
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','DGPSvDGPSISLS3sats.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','DGPSvDGPSISLS3sats.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','DGPSvDGPSISLS3sats.png'));
            end
        end

    elseif strcmpi(orbit_type,'oblate') == 1
        if SNC == 0
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/1orbitalperiod/','DGPSvDGPSISLS3sats.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/24hours/','DGPSvDGPSISLS3sats.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/10days/','DGPSvDGPSISLS3sats.png'));
            end
        elseif SNC == 1
            if ToF == 1
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/1orbitalperiod/','DGPSvDGPSISLS3sats.png'));
            elseif ToF == 24
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/24hours/','DGPSvDGPSISLS3sats.png'));
            elseif ToF == 10
                saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/snc ON/10days/','DGPSvDGPSISLS3sats.png'));
            end
        end
    end
end