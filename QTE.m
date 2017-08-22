function [EulerAngle]=QTE(q)
roll=atan2(2*(q(4)*q(1)+q(2)*q(3)),(1-2*(q(1)^2+q(2)^2))) ;
pitch=asin(2*(q(4)*q(2)-q(3)*q(1)));
yaw=atan2(2*(q(4)*q(3)+q(1)*q(2)),(1-2*(q(2)^2+q(3)^2))) ;
EulerAngle=[roll; pitch; yaw];
end