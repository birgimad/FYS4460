close all
clear all
clc

%Import data and define closed-form solution
%Data for ordered spins (up = 1) and T = 1.0
%Initial data
fileID = fopen('Datafile_velocities7_0.txt','r');
data1 = textscan( fileID, '%f %f %f') ;
fclose(fileID);
vxi = data1{1};
vyi = data1{2};
vzi = data1{3};

%final data
fileID = fopen('Datafile_velocities7_80000.txt','r');
data2 = textscan( fileID, '%f %f %f') ;
fclose(fileID);
vxf = data2{1};
vyf = data2{2};
vzf = data2{3};

%Plot histograms
% figure
% xbins = -15:1:15;
% [f,x] = hist(vxi,xbins)
% dx = diff(x(1:2));
% bar(x,f/sum(f*dx))
% 
% legend('Initial speed (x)')
% xlabel('Speed', 'fontsize',14) % x-axis label
% ylabel('Probability','fontsize',14) % y-axis label
% 
% figure
% xbinsf = -60:5:60;
% [ff,xf] = hist(vxf,xbinsf)
% dxf = diff(xf(1:2));
% bar(xf,ff/sum(ff*dxf),'r')
% 
% legend('Final speed (x)')
% xlabel('Speed', 'fontsize',14) % x-axis label
% ylabel('Probability','fontsize',14) % y-axis label

figure
xbins = -15:5:15;
[f,x] = hist(vxi,xbins)
dx = diff(x(1:2));
bar(x,f/sum(f*dx))
hold on
xbinsf = -60:5:60;
[ff,xf] = hist(vxf,xbinsf)
dxf = diff(xf(1:2));
bar(xf,ff/sum(ff*dxf),'r')

legend('Initial speed (x)','Final speed (x)')
xlabel('Speed', 'fontsize',14) % x-axis label
ylabel('Probability','fontsize',14) % y-axis label

figure
xbins = 0:2:20;
[f,x] = hist(sqrt(vxi.*vxi+vyi.*vyi+vzi.*vzi),xbins)
dx = diff(x(1:2));
bar(x,f/sum(f*dx))
hold on
xbinsf = 0:2:20;
[ff,xf] = hist(sqrt(vxf.*vxf+vyf.*vyf+vzf.*vzf),xbinsf)
dxf = diff(xf(1:2));
bar(xf,ff/sum(ff*dxf),'r')

legend('Initial speed','Final speed')
xlabel('Speed', 'fontsize',14) % x-axis label
ylabel('Probability','fontsize',14) % y-axis label



