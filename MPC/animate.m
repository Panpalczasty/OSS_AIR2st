% Shows an animation of robot state
%
% Parameters:
%   id - id of figure
%   pos - position of robot (as array)
%   ps - sim & model parameter struct
function res = animate(id, pos, ps)

    N = length(pos(:,1));
    figure(id);
    grid();
    hold on;
    plot(0,0,"k*",'LineWidth',3);
    trackline = plot(pos(1,1), pos(1,2), 'r-');
    j2plot = plot([pos(1,3) pos(1,1)],[pos(1,4) pos(1,2)], 'k*-', 'LineWidth',3);
    j1plot = plot([0 pos(1,3)],[0 pos(1,4)],'k','LineWidth',3);
    axis([-2 2 -2 2]);
    daspect([1 1 1]);

    xlabel("X [m]");
    ylabel("Y [m]");
    title("Robot Animation")
    
    x01 = ps.x0(1); 
    x02 = ps.x0(2);
    
    xf1 = ps.xf(1);
    xf2 = ps.xf(2);

    filename = sprintf("gifs/%.2f_%.2f_to_%.2f_%.2f.gif", x01, x02, xf1, xf2);

    for i = 1:1:N
        trackline.XData = pos(1:i,1);
        trackline.YData = pos(1:i,2);
        j1plot.XData = [0 pos(i,3)];
        j1plot.YData = [0 pos(i,4)];
        j2plot.XData = [pos(i,3) pos(i,1)];
        j2plot.YData = [pos(i,4) pos(i,2)];
        drawnow 
        
        frame = getframe(id);
        img = frame2im(frame);
        [img, cm] = rgb2ind(img,256);

        if i == 1
            imwrite(img, cm, filename, 'gif', 'LoopCount',inf);
        else
            imwrite(img, cm, filename, 'gif', "WriteMode","append", 'DelayTime', 1/(2*ps.sfreq));
        end
    end
    
    res = 0;
end
