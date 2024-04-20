%% This script compares two gravity solutions (usually consecutive solutions) in order to
% show the general changes between them, such as:
% - Missing stations from one solution to the next (if any)
% - stations whose gravity estimate changed by more than 0.5 mGal
% - Number of stations added
% - Changes in the mean uncertainty level of the solutions
% - Changes in the min and max uncertainty levels

% It is also possible to compare 2 solutions with the same number of
% stations. This helps to visualize changes applied to one solution w.r.t 
% its previous version. For example, effects of using one tidal model or
% another.

% Franco Sobrero, OSU, 2024

clc
clear all
close all

s1 = readtable('../gravity/gravity_Bolivia/Yr10_Franco_2020/solution_12-22-2020.txt');
s2 = readtable('../gravity/gravity_Bolivia/Yr12_Franco_2024/solution_29-Mar-2024.txt');
% s2 = readtable('../gravity/gravity_Bolivia/Yr12_Franco_2024/solution_28-Feb-2024.txt');

% s1 = readtable('../gravity/gravity_Colombia/Yr2_2022/solution_30-01-2023.txt');
% s2 = readtable('../gravity/gravity_Colombia/Yr3_2023/solution_28-Mar-2024.txt');

[~, stn2] = ismember(s1.Station, s2.Station);
[~, stn1] = ismember(s2.Station, s1.Station);

fprintf('Missing stations from S1 in S2\n')
disp(s1.Station(stn2 == 0))

r = s1.Gravity__mGal_(stn2 ~= 0) - s2.Gravity__mGal_(stn1 ~= 0);

[minr, minp] = min(r);
[maxr, maxp] = max(r);

fprintf('Min %.3f at %s max %.3f at %s\n\n', minr, s1.Station{stn1(minp)}, max(r), s1.Station{stn2(maxp)})

fprintf('Stations with abs(differences) > 0.5 mGal\n')

stn = s1.Station(stn2 ~= 0);

for i = 1:length(r)
    if abs(r(i)) > 0.5
        fprintf('%s %.3f\n\n', stn{i}, r(i))
    end
end

fprintf('Number of new stations:')
disp(length(stn1)-length(stn2))

fprintf('S1:\n')
fprintf(['Mean uncertainty: ', num2str(mean(s1.Uncertainty__mGal_)), ' mGal\n'])
fprintf(['Min uncertainty: ', num2str(min(s1.Uncertainty__mGal_)), ' mGal\n'])
fprintf(['Max uncertainty: ', num2str(max(s1.Uncertainty__mGal_)), ' mGal\n\n'])

fprintf('S2:\n')
fprintf(['Mean uncertainty: ', num2str(mean(s2.Uncertainty__mGal_)), ' mGal\n'])
fprintf(['Min uncertainty: ', num2str(min(s2.Uncertainty__mGal_)), ' mGal\n'])
fprintf(['Max uncertainty: ', num2str(max(s2.Uncertainty__mGal_)), ' mGal\n'])


%% This is just to compare 2 solutions with EXACTLY the same stations

if length(stn1) == length(stn2)
    central_latitude=mean(s1.Latitude__deg_);
    lon_bnch=s1.Longitude__deg_;
    lat_bnch=s1.Latitude__deg_;
    dot_size = 30;
    
    figure(1)
    histogram(s1.Gravity__mGal_-s2.Gravity__mGal_)
    grid on
    title('Difference between OLD and NEW Solutions')
    subtitle(['Mean: ', num2str(mean(s1.Gravity__mGal_-s2.Gravity__mGal_)), ' RMS: ', num2str(rms(s1.Gravity__mGal_-s2.Gravity__mGal_)), ' Max: ', num2str(max(abs(s1.Gravity__mGal_-s2.Gravity__mGal_)))])
    grid on
    set(gca,'FontSize',14)

    figure(2)
    subplot(1,4,1)
    plot(s1.Gravity__mGal_-s2.Gravity__mGal_)
    title('Gravity difference','FontSize',14)
    ylabel('[mGal]','FontSize',14)
    grid on
    
    subplot(1,4,2)
    plot(s1.Uncertainty__mGal_-s2.Uncertainty__mGal_)
    title('Uncertainty difference','FontSize',14)
    grid on
    
    subplot(1,4,3)
    plot(s1.Free_Air_Anomaly__mGal_-s2.Free_Air_Anomaly__mGal_)
    title('FA anom difference','FontSize',14)
    grid on
    
    subplot(1,4,4)
    plot(s1.Bouguer_Anomaly__mGal_-s2.Bouguer_Anomaly__mGal_)
    title('Bouger anom difference','FontSize',14)
    grid on
    
    figure(3)
    subplot(1,2,1)
    scatter(lon_bnch,lat_bnch,dot_size,s2.Uncertainty__mGal_,'filled')
    dar = rectangular_projection(central_latitude,'degrees');
    set(gca,'DataAspectRatio',dar)
    colormap(flipud(hot))
    Hcb=colorbar;
    colorbartitle(Hcb,'[mGal]',14)
    set(gca,'Color',0.8*[1 1 1])   % light grey background to entire plot
    set(gca,'XLim',[min(lon_bnch)-2 max(lon_bnch)+2],'YLim',[min(lat_bnch)-1 max(lat_bnch)+1])
    set(gca,'Box','on')
    hold on
    gr=0.25*[1 1 1]; % dark grey
    borders('Colombia','color',gr,'linewidth',1)
    borders('Bolivia','color',gr,'linewidth',1)
    grid on
    title('Uncertainty NEW solution')
    
    subplot(1,2,2)
    scatter(lon_bnch,lat_bnch,dot_size,s1.Uncertainty__mGal_,'filled')
    dar = rectangular_projection(central_latitude,'degrees');
    set(gca,'DataAspectRatio',dar)
    colormap(flipud(hot))
    Hcb=colorbar;
    colorbartitle(Hcb,'[mGal]',14)
    set(gca,'Color',0.8*[1 1 1])   % light grey background to entire plot
    set(gca,'XLim',[min(lon_bnch)-2 max(lon_bnch)+2],'YLim',[min(lat_bnch)-1 max(lat_bnch)+1])
    set(gca,'Box','on')
    hold on
    gr=0.25*[1 1 1]; % dark grey
    borders('Colombia','color',gr,'linewidth',1)
    borders('Bolivia','color',gr,'linewidth',1)
    grid on
    title('Uncertainty OLD solution')
end