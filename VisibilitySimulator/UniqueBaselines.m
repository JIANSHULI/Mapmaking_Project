close all; clear all;

positions = load('Compact_Omniscope.dat');


baselineMap = containers.Map();
for i = 1:length(positions)
    for j = i:length(positions)
        if i ~= j
            deltaSouth = -(positions(i,2) - positions(j,2));
            deltaEast = positions(i,1) - positions(j,1);
            key = [num2str(deltaSouth) ',' num2str(deltaEast)];
            if isKey(baselineMap,key)
                baselineMap(key) = baselineMap(key) + 1;
            else
                baselineMap(key) = 1;
            end
        end
    end
end

positions = positions/2;

% for i = 1:length(positions)
%     for j = i:length(positions)
%         if i ~= j
%             deltaSouth = -(positions(i,2) - positions(j,2));
%             deltaEast = positions(i,1) - positions(j,1);
%             key = [num2str(deltaSouth) ',' num2str(deltaEast)];
%             if isKey(baselineMap,key)
%                 baselineMap(key) = baselineMap(key) + 1;
%             else
%                 baselineMap(key) = 1;
%             end
%         end
%     end
% end

baselineKeys = keys(baselineMap);
baselineValues = values(baselineMap);

for k = 1:length(baselineKeys)
    commaPos = strfind(baselineKeys{k},',');
    deltaSouth= str2num(baselineKeys{k}(1:(commaPos-1)));
    deltaEast = str2num(baselineKeys{k}((commaPos+1):end));
    baselineMatrix(k,:) = [ deltaSouth deltaEast baselineValues{k}];
end

baselineMatrix = flipud(sortrows(baselineMatrix,3));

filename = 'IntendedUniqueBaselines.dat';
fileID = fopen(filename,'w');
for (n = 1:size(baselineMatrix,1))
    fprintf(fileID,'%f %f %f %d',baselineMatrix(n,1),baselineMatrix(n,2),0,baselineMatrix(n,3));
    if n < size(baselineMatrix,1)
        fprintf(fileID,'\n')
    end
end
fclose(fileID);