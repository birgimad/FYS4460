close all
clear all
clc

%Import data and define closed-form solution
%Data for ordered spins (up = 1) and T = 1.0
filename = 'Initial_velocity_gaussian_dist.xlsx';
sheet = 2;
xlRange = 'A3:C2050';

[v,T,vT] = xlsread(filename, sheet, xlRange);
vx=v(:,1);
vy=v(:,2);
vz=v(:,3);

y = [-600:30:600];
norm = normpdf(y,0,200);

%Plot histograms
figure
xbins = -600:30:600;
[f,x] = hist(vy,xbins)
dx = diff(x(1:2));
bar(x,f/sum(f*dx))
hold on
plot(y,norm,'r','LineWidth',2)
legend('Generated velocity, v_y','Gaussian dist.')
xlabel('Velocity', 'fontsize',14) % x-axis label
ylabel('Probability','fontsize',14) % y-axis label




