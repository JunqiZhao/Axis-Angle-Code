function [q] = invQuat(q)
q = [-1*q(1:3);q(4)];
end
