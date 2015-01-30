clear all; 
close all;

MWA = load('mwa_dataX');
for f = 1:9
    data = reshape(MWA(((1+(f-1)*45*180):(f*45*180)),7),45,180);
    dlmwrite(['mwa_beam_xx_' num2str(100+10*f) '.dat'],data,' ');
end

MWA = load('mwa_dataY');
for f = 1:9
    data = reshape(MWA(((1+(f-1)*45*180):(f*45*180)),7),45,180);
    dlmwrite(['mwa_beam_yy_' num2str(100+10*f) '.dat'],data,' ');
end
