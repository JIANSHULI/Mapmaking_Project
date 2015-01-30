close all; clear all;
freqs = {'125.586', '150.195', '175.196'};
angularResolution = '2';


for n = 1:length(freqs)
    file = dir(['Results/*' freqs{n} '*' angularResolution '_Deg_*']);
    facet = load(['Results/' file.name]);
    figure(n); imagesc(log10(facet)); colorbar
    title([freqs{n} ' MHz'])
    xlabel('RA bin'); ylabel('Dec bin');
end
