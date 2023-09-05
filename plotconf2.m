       
function    [output] = plotconf2(Xhcw)
            %formation 2
            origin = plot3(0,0,0,'kd'); %virtual chief placed in the origin.
            sat1 = plot3(Xhcw(1,:),Xhcw(2,:),Xhcw(3,:),'b-','LineWidth',1); %trajectory of sat1, if not using a gif.
            sat2 = plot3(Xhcw(7,:),Xhcw(8,:),Xhcw(9,:),'m-','LineWidth',1); %trajectory of sat2, if not using a gif.
            sat3 = plot3(Xhcw(13,:),Xhcw(14,:),Xhcw(15,:),'r-','LineWidth',1); %trajectory of sat3, if not using a gif.
            sat4 = plot3(Xhcw(19,:),Xhcw(20,:),Xhcw(21,:),'b-','LineWidth',1); %trajectory of sat1, if not using a gif.

            %initial and final positions
            initial_pos1 = plot3(Xhcw(1,1), Xhcw(2,1), Xhcw(3,1), 'sr', 'MarkerFaceColor', 'k','MarkerSize',10); %deputy initial position
            final_pos1   = plot3(Xhcw(1,end), Xhcw(2,end), Xhcw(3,end), 'sr', 'MarkerFaceColor', 'r'); %formation final position
            ite_sat1 = scatter3(Xhcw(3,1), Xhcw(2,1), Xhcw(3,1),[20],'or','filled','MarkerEdgeColor','k');  %first iteration for our formation

            initial_pos2 = plot3(Xhcw(7,1), Xhcw(8,1), Xhcw(9,1), 'sb', 'MarkerFaceColor', 'k','MarkerSize',10); %deputy initial position
            final_pos2   = plot3(Xhcw(7,end), Xhcw(8,end), Xhcw(9,end), 'sb', 'MarkerFaceColor', 'r'); %formation final position
            ite_sat2 = scatter3(Xhcw(7,1), Xhcw(8,1), Xhcw(9,1),[20],'ob','filled','MarkerEdgeColor','k');  %first iteration for our formationinitial_pos = plot3(Xhcw(1,1), Xhcw(7,1), Xhcw(13,1), 'sr', 'MarkerFaceColor', 'b','MarkerSize',10); %deputy initial position

            initial_pos3 = plot3(Xhcw(13,1), Xhcw(14,1), Xhcw(15,1), 'sg', 'MarkerFaceColor', 'k','MarkerSize',10); %deputy initial position
            final_pos3   = plot3(Xhcw(13,end), Xhcw(14,end), Xhcw(15,end), 'sg', 'MarkerFaceColor', 'r'); %formation final position
            ite_sat3 = scatter3(Xhcw(13,1), Xhcw(14,1), Xhcw(15,1),[20],'og','filled','MarkerEdgeColor','k');  %first iteration for our formation

            initial_pos4 = plot3(Xhcw(19,1), Xhcw(20,1), Xhcw(21,1), 'sm', 'MarkerFaceColor', 'k','MarkerSize',10); %deputy initial position
            final_pos4  = plot3(Xhcw(19,end), Xhcw(20,end), Xhcw(21,end), 'sm', 'MarkerFaceColor', 'r'); %formation final position
            ite_sat4 = scatter3(Xhcw(19,1), Xhcw(20,1), Xhcw(21,1),[20],'om','filled','MarkerEdgeColor','k');  %first iteration for our formation

            output = [origin sat1 sat2 sat3 sat4 initial_pos1 final_pos1 ite_sat1 initial_pos2 final_pos2 ite_sat2 initial_pos3 final_pos3 ite_sat3 initial_pos4 final_pos4 ite_sat4];