function [ mu1, sigma1, w ] = myekf( mu0, sigma0, u1, zz, errZ, sigr, sigw, dt, i)

%% Initial Setup
w0 = u1(i,:)'; %turning rate measurement at time0
w1 = u1(i+1,:)'; %turning rate measurement at time1
q_hat_avg_0 = mu0(1:4)'; %initial quaternion
b0 = mu0(5:7); %initial bias
w_hat_0 = w0 - b0; %initial turning rate measurement with corrected bias
%% State Prediction
b1 = b0; %assuming bias constant
w_hat_1 = w1 - b1; %estimate turning rate at next timepoint
w_avg = (w_hat_0 + w_hat_1) / 2;
q_hat_avg_1 = NormalizeV((exp(1/2*Omega(w_avg)*dt) + 1/48*(Omega(w_hat_1)*Omega(w_hat_0) - Omega(w_hat_0)*Omega(w_hat_1))*dt^2) * q_hat_avg_0); % estimating quaternion in next time point

%% Error Propagation during State Prediction
% Computation of transition matrix-Phi
if norm(w_avg) < 0.00001 %Too small, may induce numerical instability. Estimating instead.
    Theta = eye(3) - dt*skew(w_avg) + ((dt.^2)/2)*(skew(w_avg)*skew(w_avg));
    Psi = -eye(3)*dt + ((dt^2)/2)*skew(w_avg) - ((dt^3)/6)*(skew(w_avg)*skew(w_avg));
else
    Theta = cos(norm(w_avg)*dt)*eye(3) - sin(norm(w_avg)*dt)*skew(w_avg/norm(w_avg)) + (1 -cos(norm(w_avg)*dt))*(w_avg/norm(w_avg))*(w_avg'/norm(w_avg));
    Psi = -eye(3)*dt + (1/norm(w_avg).^2)*(1-cos(norm(w_avg)*dt))*skew(w_avg) - (1/norm(w_avg).^3)*(norm(w_avg)*dt - sin(norm(w_avg)*dt))*(skew(w_avg)*skew(w_avg));
end
Phi = [Theta Psi; zeros(3) eye(3)];
% Computation of process noise Qd
if norm(w_avg) < 0.00001  %Too small, may induce numerical instability. Estimating instead.
    Q11 = (sigr.^2)*dt*eye(3) + (sigw.^2)*(eye(3)*dt^3/3 + (dt^5/60)*(skew(w_avg)*skew(w_avg)));
    Q12 = -(sigw^2) * ( eye(3)*dt.^2/2 - (dt^3/6)*skew(w_avg) + (dt^4/24)*(skew(w_avg)*skew(w_avg)));
else
    Q11 = (sigr.^2)*dt*eye(3) + (sigw.^2)*( eye(3)*dt.^3/3 + (((norm(w_avg)*dt)^3/3 + 2*sin(norm(w_avg)*dt) - 2*norm(w_avg)*dt )/ (norm(w_avg)^5))*(skew(w_avg)*skew(w_avg)));
    Q12 = -(sigw^2) * ( eye(3)*dt.^2/2 - ((norm(w_avg)*dt - sin(norm(w_avg)*dt))/(norm(w_avg)^3))*skew(w_avg) + (((norm(w_avg)*dt)^2/2 + cos(norm(w_avg)*dt) - 1)/(norm(w_avg)^4))*(skew(w_avg)*skew(w_avg)));
end
Q22 = (sigw^2)*dt*eye(3);
Qd = [Q11 Q12; Q12' Q22];
% Compute the state covariance matrix according to the Extended Kalman Filter equation
sigma1_ = Phi*sigma0*Phi' + Qd;

%% Computation of the Kalman gain
z_hat     = convertToAxisAngle(q_hat_avg_1);
% Compute the measurement matrix H
H = [skew(q_hat_avg_1) zeros(3)];
% Compute residual r according to r = z - z_hat, while transforming the measurement into quaternion
z_quat = convertToQuaternion(zz);
z_hat_quat = convertToQuaternion(z_hat);
r = multiplyQuaternion(z_quat,invQuat(z_hat_quat));
% Compute the covariance of the residual S as
S = H*sigma1_*H' + errZ;
% Compute the Kalman gain K
K = sigma1_*(H'*inv(S));

%% Measurement Update - State
% Compute the correction deltaX = [2*dq;db]
deltaX = K*r(1:3);
dq = deltaX(1:3)/2;
db = deltaX(4:6);
% Update the quaternion
if (dq'*dq) > 1
    dq_hat_avg_1 = (1/sqrt(1 + (dq'*dq))) * [dq ; 1];
else
    dq_hat_avg_1 = [dq ; (sqrt(1 - (dq'*dq)))];
end
q = multiplyQuaternion(dq_hat_avg_1, q_hat_avg_1);
% Update the bias
b = b1' + db;
mu1 = [q;b]';
% Update the estimated turn rate using the new estimate for the bias
w = w1 - b;

%% Measurement Update - Error
% Compute the new updated Covariance matrix
sigma1 = (eye(6) - K*H) * sigma1_ * (eye(6) - K*H)' + K*errZ*K';
end