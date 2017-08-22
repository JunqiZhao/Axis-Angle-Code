%% EXTENDED KALMAN FILTER WITH QUATERNIONS
clear all
close all
clc

%% Initial Setup outsite the EKF
mydata=xlsread('Roll Test.xlsx');
u1=mydata(:,5:7);%gyro measurement
mu0=(mydata(1,12:18)); %initial state vector
dt=0.1; %delta T
sigr  = 0.01; % rate noise covariance
sigw  = 0.01; % gyro bias covariance
sigma0=eye(6)*10^-4; % initial state covariance matrix
G_x = mydata(:,8);
G_y = mydata(:,9);
G_z = mydata(:,10); 
%% Buffer 
history.mu1    = zeros(length(mydata),7);
history.sigma1 = zeros(6,6);
history.w      = zeros(length(mydata),3);
history.Euler  = zeros(3,length(mydata));
history.EulerZ  = zeros(3,length(mydata));
history.Axis   = zeros(4,length(mydata));
fai_acc        = zeros(1,length(mydata));

z=zeros(4,length(mydata));
%% EFK Execution
for i =1:length(mydata)-1
   %% Assuming using the Axis angle
   z(:,i)=[mydata(i,5)/sqrt(mydata(i,5)^2+mydata(i,6)^2+mydata(i,7)^2),mydata(i,6)/sqrt(mydata(i,5)^2+mydata(i,6)^2+mydata(i,7)^2),mydata(i,7)/sqrt(mydata(i,5)^2+mydata(i,6)^2+mydata(i,7)^2),dt*sqrt(mydata(i,5)^2+mydata(i,6)^2+mydata(i,7)^2)]';
   zz=z(4,i);
   errZ=[2.303667158010755e-05, 0, 0; 0, 2.303667158010755e-05, 0; 0, 0, 2.303667158010755e-05]; % This is measurement noise, and related to the measurement method, here we assuming using the IMU sensor and used the noise matrix from IMU
   [mu1, sigma1, w ] = myekf( mu0, sigma0, u1, zz, errZ, sigr, sigw, dt, i);
   %% Store the result
   history.mu1(i,:)  = mu1';
   history.w(i,:)    = w';
   history.Euler(:,i)= QTE(history.mu1(i,1:4));
   history.EulerZ(:,i)= QTE(convertToQuaternion(zz));
   history.fai_acc(1,i) = atan2(-G_y(i),G_z(i));
   sigma0=sigma1;
end
plot(history.fai_acc(1,:)) 
hold on 
plot(history.Euler(1,:))
hold off

