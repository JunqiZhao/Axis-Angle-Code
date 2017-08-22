%% EXTENDED KALMAN FILTER WITH QUATERNIONS
clear all
close all
clc

%% Initial Setup outsite the EKF
mydata=xlsread('Test0808.xlsx')
u1=mydata(:,5:7);%initial measurement
mu0=(mydata(1,12:18)); %initial state vector
dt=0.01 %delta T
sigr  = 0.1; % rate noise covariance (measurement noise)
sigw  = 0.1; % gyro bias covariance
sigma0=zeros(6,6) % initial state covariance matrix
%% Buffer 
history.mu1    = zeros(length(mydata),7);
history.sigma1 = zeros(6,6);
history.w      = zeros(length(mydata),3);
%% EFK Execution
for i=1:length(mydata)
   %% Assuming using the Axis angle
   z= z=[1,0,0,mydata(i,5)]';
   [mu1, sigma1, w ] = myekf( mu0, sigma0, u1, z, errZ, sigr, sigw, dt)
   %% Store the result
   history.mu1(i,:) = mu1';
   history.w(i,:)   = w';
   mu0=mu1
   sigma0=sigma1
end
