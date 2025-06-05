% SAVEFIGS example

%% Script
% Cleaning
close all;

% Plot something
x = linspace(0,2*pi,100);
y1 = sin(x);
y2 = tan(x);
figure(2); plot(x,y1);
figure(3); plot(x,y2);

% I'm choosing a save folder (if I don't they will be saved into .\fig\
savedir = 'myfigures';

% I can just save the figures with their numbers, so that the created name
% would be "fig2", "fig3"...
saveFigs(savedir);

% or I can specify the name of the figures
f2 = figure(2);
f2.Name = 'Sine plot';
f3 = figure(3);
f3.Name = 'Tan plot';

saveFigs(savedir);

% I can choose the format of the output files as a Param/Value option, and
% also one particular figure to save
saveFigs(savedir,'format','pdf','handle',f3);