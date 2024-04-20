clc
clear all

data = readtable('../data/data_Bolivia/5_September_2023_Uyuni/Coordenadas_GPS_Eric/Eric_Coords_with_height_correction/20231201_Coords_update_Uyuni2023_wo_tie_stations');

load('../w1_bolivia/database.mat');
updated = 0;
not_found = 0;
conflict = 0;
for i = 1:size(data,1)
    
    benchmark = lower(data.Var1(i));
    
    index = find(ismember({database.benchmarks.name}, benchmark));
    
    lat = data.Var6(i);
    lon = data.Var7(i);
    
    if lat == 0 && lon == 0
        lla = ecef2lla([data.Var3(i) data.Var4(i) data.Var5(i)]);
        lat = lla(1);
        lon = lla(2);
    end
        
    if index ~= 0
        alat = database.benchmarks(index).lat;
        alon = database.benchmarks(index).lon;
        
        x = database.benchmarks(index).x;
        y = database.benchmarks(index).y;
        z = database.benchmarks(index).z;
        
        distance = 2*asin(sqrt(sin((deg2rad(alat)-deg2rad(lat))/2)^2 + cos(deg2rad(lat)) * cos(deg2rad(alat)) * sin((deg2rad(alon)-deg2rad(lon))/2)^2))*6371000;
        
        if abs(database.benchmarks(index).height) > 1e-2
            distance = sqrt((x-data.Var3(i)).^2 + (y-data.Var4(i)).^2 + (z-data.Var5(i)).^2);
            
            fprintf('WARNING! %s(%i) already has an assigned height! Difference between coordinates is %.1f km\n', benchmark{1}, i, distance./1000)
            conflict = conflict + 1;
        else
            % check apriori with final coordinate
            if distance > 2000
                fprintf('WARNING! %s(%i) Difference between apriori coordinate and final is %.1f km\n', benchmark{1}, i, distance./1000)
                conflict = conflict + 1;
            else
                bench = CssBenchmark(benchmark{1}, data.Var3(i), data.Var4(i), data.Var5(i), 0);
                database.benchmarks(index).UpdateCoords(bench);
                updated = updated + 1;
            end
        end
    else
        not_found = not_found + 1;
        
        % not found, find a candidate using distance
        alat = [database.benchmarks.lat]';
        alon = [database.benchmarks.lon]';
        
        distance = 2*asin(sqrt(sin((deg2rad(alat)-deg2rad(lat))./2).^2 + cos(deg2rad(lat)) .* cos(deg2rad(alat)) .* sin((deg2rad(alon)-deg2rad(lon))./2).^2)).*6371000;
        
        [md, k] = min(distance);
        
        fprintf(' -- Could not find %s(%i) in the database: closest (%.3f km) is %s\n', benchmark{1}, i, md/1000, database.benchmarks(i).name)
    end
    
end

j = 1;
zero_h = {};

for i = 1:length(database.benchmarks)
    
    if abs(database.benchmarks(i).height) < 1e-2 && isempty(database.benchmarks(i).absolute_g)
        zero_h{j,1} = database.benchmarks(i).name;
        j = j+1;
    end
end

fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('Done updating coordinates: %i updated, %i not found, %i with conflicts\n', updated, not_found, conflict)
fprintf('List of benchmarks with zero height saved to zero_h variable\n')
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')


