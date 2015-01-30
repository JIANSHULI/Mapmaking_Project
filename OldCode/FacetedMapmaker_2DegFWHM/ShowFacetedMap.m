close all; clear all;

map = load('map.dat');
figure(1);
imagesc(map'); colorbar;
set(gca,'YDir','normal')
xlabel('RA Bin')
ylabel('Dec Bin')

mapFull = load('unnormalized_full_map.dat');
figure(2);
imagesc(mapFull'); colorbar;
set(gca,'YDir','normal')
xlabel('RA Bin')
ylabel('Dec Bin')

PSF = load('PSF.dat');
figure(3);
imagesc(reshape(PSF(200,:),length(map),[]))

%%
% 
% snapshotMap = load('snapshotMap.dat');
% figure(3);
% imagesc(snapshotMap'); colorbar;
% set(gca,'YDir','normal')
% xlabel('RA Bin')
% ylabel('Dec Bin')

% %%
% facetGrid = load('facetGridCopy_Real.dat') + sqrt(-1)*load('facetGridCopy_Imag.dat');
% figure(4)
% imagesc(abs(facetGrid))
% 
% %% 
% KAt = load('KAtranspose_Real.dat') + sqrt(-1)*load('KAtranspose_Imag.dat');
% %%
% figure(5)
% imagesc(abs(KAt))
% 
% %%
% nCropped = sqrt(size(KAt,2));
% 
% croppedSnapshotMap = snapshotMap(size(snapshotMap,1)/2 + (-nCropped/2 + 1:nCropped/2), size(snapshotMap,2)/2 + (-nCropped/2 + 1:nCropped/2));
% figure(6);imagesc(croppedSnapshotMap)
% 
% 
% 
% Ninv = load('Ninv.dat');
% %figure(10); imagesc(Ninv)
% reconstructedMap = real(reshape((KAt.'*reshape(Ninv.*facetGrid,[],1)),nCropped,nCropped));
% %reconstructedMap = real(reshape((KAt.'*reshape(facetGrid,[],1)),nCropped,nCropped));
% 
% 
% 
% figure(7);imagesc(reconstructedMap)
% 
% figure(8); imagesc(croppedSnapshotMap./reconstructedMap); colorbar



