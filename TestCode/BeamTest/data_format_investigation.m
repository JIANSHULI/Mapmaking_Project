clear all; 
close all;

name = 'mwa_dataX';

MWA = load(name);

%%

for f = 1:9
    data = reshape(MWA(((1+(f-1)*45*180):(f*45*180)),7),45,180);
    %write(array2table(data), ['mwa_beamXX_' num2str(100+10*f) '.dat'], 'Delimiter', ' ');
    dlmwrite(['mwa_beam_xx_' num2str(100+10*f) '.dat'],data,' ');
end

%%

mesh(0:(360/179):360, 0:(90/44):90, data)