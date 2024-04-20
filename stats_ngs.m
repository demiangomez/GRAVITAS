clc
clear all

% stats of measured lines
%ff = dir('../w1_bolivia/lines/*.mat');
ff = dir('../w1_colombia/lines/*.mat');

array = [];

for i = 1:length(ff)
    load([ff(i).folder '/' ff(i).name])
    
    for j = 1:length(line_struct.observations)
        o = line_struct.observations(j);
        array(end+1) = decyear(o.timestamp);
    end
end

figure(1)
histogram(array, 200)
grid on
ylabel('Count')
xlabel('Epochs')
title('Number of observations')
set(gca,'FontSize',15)