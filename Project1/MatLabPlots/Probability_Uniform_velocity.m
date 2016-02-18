close all
clear all
clc

%Import data and define closed-form solution
%Data for ordered spins (up = 1) and T = 1.0
filename = 'Initial_velocity_gaussian_dist.xlsx';
sheet = 5;
xlRange = 'A1:C2048';

[v,T,vT] = xlsread(filename, sheet, xlRange);
vxi=v(:,1);
vyi=v(:,2);
vzi=v(:,3);

filename = 'Initial_velocity_gaussian_dist.xlsx';
sheet = 6;
xlRange = 'A1:C2048';

[v,T,vT] = xlsread(filename, sheet, xlRange);
vxf=v(:,1);
vyf=v(:,2);
vzf=v(:,3);

%Plot histograms
figure
xbins = -12:1:12;
[f,x] = hist(vxi,xbins)
dx = diff(x(1:2));
bar(x,f/sum(f*dx))

legend('Initial speed (x)')
xlabel('Speed', 'fontsize',14) % x-axis label
ylabel('Probability','fontsize',14) % y-axis label

figure
xbinsf = -20:1:20;
[ff,xf] = hist(vxf,xbinsf)
dxf = diff(xf(1:2));
bar(xf,ff/sum(ff*dxf),'r')

legend('Final speed (x)')
xlabel('Speed', 'fontsize',14) % x-axis label
ylabel('Probability','fontsize',14) % y-axis label

figure
xbins = -12:1:12;
[f,x] = hist(vxi,xbins)
dx = diff(x(1:2));
bar(x,f/sum(f*dx))
hold on
xbinsf = -20:1:20;
[ff,xf] = hist(vxf,xbinsf)
dxf = diff(xf(1:2));
bar(xf,ff/sum(ff*dxf),'r')

legend('Initial speed (x)','Final speed (x)')
xlabel('Speed', 'fontsize',14) % x-axis label
ylabel('Probability','fontsize',14) % y-axis label

figure
xbins = 0:1:20;
[f,x] = hist(sqrt(vxi.*vxi+vyi.*vyi+vzi.*vzi),xbins)
dx = diff(x(1:2));
bar(x,f/sum(f*dx))
hold on
xbinsf = 0:1:30;
[ff,xf] = hist(sqrt(vxf.*vxf+vyf.*vyf+vzf.*vzf),xbinsf)
dxf = diff(xf(1:2));
bar(xf,ff/sum(ff*dxf),'r')

legend('Initial speed','Final speed')
xlabel('Speed', 'fontsize',14) % x-axis label
ylabel('Probability','fontsize',14) % y-axis label



