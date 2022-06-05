function [] = saveToFile(id, folder, type, ps)
%SAVETOFILE Summary of this function goes here
%   Detailed explanation goes here
    x01 = ps.x0(1); 
    x02 = ps.x0(2);
    xf1 = ps.xf(1);
    xf2 = ps.xf(2);

    foldername = sprintf("%s/fig%.2f_%.2f_to_%.2f_%.2f", folder, x01, x02, xf1, xf2);
    if not(isfolder(foldername))
        mkdir(foldername)
    end
    
    filename = sprintf("%s/%s.png", foldername, type);
    saveas(id,filename)
end

