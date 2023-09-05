function savingplot(configuration,orbit_type,ToF,output)

global save_plot

switch configuration
    case 1
        title('Closed Orbit via HCW equations (config. 1)');
        legend([output(1:7)],'Virtual Target','Chaser1 Uncontrolled trajectory',...
            'Chaser2 Uncontrolled trajectory','Chaser3 Uncontrolled trajectory','Chasers1 Starting positions','Chasers2 Starting positions','Chasers3 Starting positions','FontSize',22)
        hold off;

        if save_plot == 1
            if strcmpi(orbit_type,'circular') == 1
                if ToF == 1 %1 hour
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/1orbitalperiod/', 'HCWconfig1.png'));
                    % exportgraphics(gcf,'/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/config1/1orbitalperiod/HCWconfig1.png');
                elseif ToF == 24 %24 hours
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/24hours/','HCWconfig1.png'));
                elseif ToF == 10 %10 days
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config1/10days/','HCWconfig1.png'));
                end
            elseif strcmpi(orbit_type,'oblate') == 1
                if ToF == 1 %1 hour
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/1orbitalperiod/', 'HCWconfig1.png'));
                    % exportgraphics(gcf,'/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/config1/1orbitalperiod/HCWconfig1.png');
                elseif ToF == 24 %24 hours
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/24hours/','HCWconfig1.png'));
                elseif ToF == 10 %10 days
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config1/10days/','HCWconfig1.png'));
                end
            end
        end

    case 2
        title('Closed Orbit via HCW equations (config. 2)');
        legend([output(1:9)],'Virtual Target','Chaser1 Uncontrolled trajectory',...
            'Chaser2 Uncontrolled trajectory','Chaser3 Uncontrolled trajectory','Chaser4 Uncontrolled trajectory','Chasers1 Starting positions','Chasers2 Starting positions','Chasers3 Starting positions',...
            'Chasers4 Starting positions','FontSize',22)
        hold off;
        if save_plot == 1
            if strcmpi(orbit_type,'circular') == 1
                if ToF == 1 %1 hour
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/1orbitalperiod/','HCWconfig2.png'));
                elseif ToF == 24 %24 hours
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/24hours/','HCWconfig2.png'));
                elseif ToF == 10 %10 days
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config2/10days/','HCWconfig2.png'));
                end
            elseif strcmpi(orbit_type,'oblate') == 1
                if ToF == 1 %1 hour
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config2/1orbitalperiod/','HCWconfig2.png'));
                elseif ToF == 24 %24 hours
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config2/24hours/','HCWconfig2.png'));
                elseif ToF == 10 %10 days
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config2/10days/','HCWconfig2.png'));
                end
            end
        end

    case 3
        title('Closed Orbit via HCW equations (config. 3)');
        legend([output(1:9)],'Chaser1 Uncontrolled trajectory',...
            'Chaser2 Uncontrolled trajectory','Chaser3 Uncontrolled trajectory','Chaser4 Uncontrolled trajectory','Central Chief','Chaser1 Starting positions','Chaser2 Starting positions','Chaser3 Starting positions',...
            'Chaser4 Starting positions','FontSize',22)
        %hold off;
        if save_plot == 1
            if strcmpi(orbit_type,'circular') == 1
                if ToF == 1 %1 hour
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/1orbitalperiod/', 'HCWconfig3.png'));
                elseif ToF == 24 %24 hours
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/24hours/', 'HCWconfig3.png'));
                elseif ToF == 10 %10 days
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/circular orbit/config3/10days/', 'HCWconfig3.png'));
                end
            elseif strcmpi(orbit_type,'oblate') == 1
                if ToF == 1 %1 hour
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config3/1orbitalperiod/', 'HCWconfig3.png'));
                elseif ToF == 24 %24 hours
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config3/24hours/', 'HCWconfig3.png'));
                elseif ToF == 10 %10 days
                    saveas(ancestor(gcf,'figure'), fullfile('/Users/ludovicogregori/Desktop/SURREY/FINAL PROJECT/matlab codes/orbit&datafilter/figures/oblate orbit/config3/10days/', 'HCWconfig3.png'));
                end
            end
        end

end
