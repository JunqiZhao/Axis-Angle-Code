
w_avg=[1,2,3,1]';
dt=0.01;
w_hat_1=[2,1,3,1]';
w_hat_0=[1,2,2,2]';
q_hat_avg_0=[1,0,0,0]';
q_hat_avg_1 = NormalizeV((exp(1/2*Omega(w_avg)*dt) + 1/48*(Omega(w_hat_1)*Omega(w_hat_0) - Omega(w_hat_0)*Omega(w_hat_1))*dt^2) * q_hat_avg_0);
