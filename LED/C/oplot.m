% octave file to plot LED data
clear all; close all; clc;
%% read in data
U=dlmread('LEDOutput.csv',' ');
p=dlmread('LEDInput.txt',' ',[1,0,1,4])
x=p(1)
dx=p(2)
t=p(3)
dt=p(4)
a=p(5)

figure(1)
[m,n]=size(U)
[X,T]=meshgrid(0:dx:x,0:dt:t);
surf(X,T,U(:,3:n-3),'LineWidth',0)
xlabel('Length (x)}','interpreter','tex')
ylabel('time (t)'),zlabel('u(x,t) value');
