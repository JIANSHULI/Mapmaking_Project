close all; clear all;

%2 * 1.3806488e-23 * (150.195*1e6)^2 / (299792458)^2 * 1e26
tempFactor = 693.079031772476355399703606963;
n = 1;
prefixes = {''};

trueSky = load([prefixes{n} 'trueSky.dat']);
PSF = load([prefixes{n} 'PSF.dat']);
map = load([prefixes{n} 'map.dat']);
coords = load([prefixes{n} 'pixelCoordinates.dat']) * 360/2/pi;
coordsExtended = load([prefixes{n} 'extendedPixelCoordinates.dat']) * 360/2/pi;
Dmatrix = load([prefixes{n} 'Dmatrix.dat']);

ras = coords(:,1);
decs = coords(:,2);
rasExtended = coordsExtended(:,1);
decsExtended = coordsExtended(:,2);
nPixels = length(ras);
nPixelsExtended = length(rasExtended);
centralPSF = PSF(round((nPixels+1)/2),:);
convolvedSky = PSF*trueSky;


%trueSky = zeros(size(trueSky));
%extendedPixels = load([prefixes{n} 'extendedHealpixPixels.dat']);
%trueSky = (extendedPixels == 145530);
%trueSky = (extendedPixels == 8342);


%%
mSize = 100;

figure(n); clf
set(n,'position',[ 221         224        1218         650])
ha = tight_subplot(2,3,[.1 .05],[.05 .05],[.05 .05]);


axes(ha(1)); 
scatter(rasExtended, decsExtended, mSize, trueSky, 'fill','Marker','d','MarkerEdgeColor','none');
axis square;
colorbar; title('True Sky');
xlabel('Right Ascension'); ylabel('Declination');
mapRange = get(gca,'CLim');
raRange = get(gca,'XLim');
decRange = get(gca,'YLim');

%catalog = load([prefixes{n} 'pointSources.dat']);
%hold on; scatter(catalog(:,1) * 360/2/pi, catalog(:,2) * 360/2/pi, 20, catalog(:,5), 'fill','Marker','o','MarkerEdgeColor','none'); hold off;
%caxis([0 5000])

axes(ha(2));
scatter(ras, decs, mSize, map, 'fill','Marker','d','MarkerEdgeColor','none');
colorbar; title('Dirty Map');
xlabel('Right Ascension'); ylabel('Declination');
mapRange = get(gca,'CLim');
axis square
set(gca,'XLim',raRange); set(gca,'YLim',decRange);

axes(ha(3));
scatter(rasExtended, decsExtended, mSize, centralPSF, 'fill','Marker','d','MarkerEdgeColor','none');
colorbar; title('Central PSF');
xlabel('Right Ascension'); ylabel('Declination');
mapRange = get(gca,'CLim');
axis square
set(gca,'XLim',raRange); set(gca,'YLim',decRange);

axes(ha(4));
scatter(ras, decs, mSize, convolvedSky, 'fill','Marker','d','MarkerEdgeColor','none');
colorbar; title('Convoled True Sky');
xlabel('Right Ascension'); ylabel('Declination');
mapRange = get(gca,'CLim');
axis square
set(gca,'XLim',raRange); set(gca,'YLim',decRange);

axes(ha(5));
scatter(ras, decs, mSize, map-convolvedSky, 'fill','Marker','d','MarkerEdgeColor','none');
colorbar; title('Dirty Map - Convolved Sky');
xlabel('Right Ascension'); ylabel('Declination');
mapRange = get(gca,'CLim');
axis square
set(gca,'XLim',raRange); set(gca,'YLim',decRange);

axes(ha(6));
scatter(ras, decs, mSize, map-convolvedSky, 'fill','Marker','d','MarkerEdgeColor','none');
colorbar; title('Dirty Map - Convolved Sky');
xlabel('Right Ascension'); ylabel('Declination');
mapRange = get(gca,'CLim');
axis square
set(gca,'XLim',raRange); set(gca,'YLim',decRange);




mapRange1 = get(ha(2),'CLim');
mapRange2 = get(ha(4),'CLim');
mapRange3 = get(ha(5),'CLim');
%mapRange = [min([mapRange1(1),mapRange2(1),mapRange3(1)]), max([mapRange1(2),mapRange2(2),mapRange3(2)])];
mapRange = mapRange2;
%axes(ha(1)); set(gca,'CLim',mapRange);

%axes(ha(2)); set(gca,'CLim',mapRange);
axes(ha(4)); set(gca,'CLim',mapRange);
axes(ha(5)); set(gca,'CLim',mapRange);



%set(n,'Color',[1 1 1])
%export_fig(gcf,sprintf('Angular_Resolution_%f_Degrees.png',angularResolutions(n)),'-nocrop','-r200')


disp(['Map vs. Convolved Sky Error is ' num2str(norm(map-convolvedSky)/norm(convolvedSky))])
