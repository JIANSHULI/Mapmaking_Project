close all; clear all;

%2 * 1.3806488e-23 * (150.195*1e6)^2 / (299792458)^2 * 1e26
tempFactor = 693.079031772476355399703606963;
n = 1;
prefixes = {'/Volumes/EoR/Storage/Mapmaking_Results/Results_Pictures/'};

trueSky = load([prefixes{n} 'trueSky.dat']);

map = load([prefixes{n} 'map.dat']);
coords = load([prefixes{n} 'pixelCoordinates.dat']) * 360/2/pi;
coordsExtended = load([prefixes{n} 'extendedPixelCoordinates.dat']) * 360/2/pi;

ras = coords(:,1);
decs = coords(:,2);
rasExtended = coordsExtended(:,1);
decsExtended = coordsExtended(:,2);
nPixels = length(ras);
nPixelsExtended = length(rasExtended);


%%
mSize = 100;

figure(n); clf
set(n,'position',[ 221         224        1218         650])
ha = tight_subplot(1,2,[.1 .05],[.05 .05],[.05 .05]);


axes(ha(1)); 
scatter(rasExtended, decsExtended, mSize, trueSky, 'fill','Marker','d','MarkerEdgeColor','none');
axis square;
colorbar; title('True Sky');
xlabel('Right Ascension'); ylabel('Declination');
mapRange = get(gca,'CLim');
raRange = get(gca,'XLim');
decRange = get(gca,'YLim');
 set(gca,'Color',[0 0 0])
    set(gca,'XColor',[1 1 1],'YColor',[1 1 1])

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
 set(gca,'Color',[0 0 0])
    set(gca,'XColor',[1 1 1],'YColor',[1 1 1])

colormap hot
set(gcf,'Color',[0 0 0])