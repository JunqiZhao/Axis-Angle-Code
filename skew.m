function [ qx ] = skew( q )   
    qx = [ 0   -q(3)  q(2);
          q(3)   0   -q(1);
         -q(2)  q(1)   0 ];
end
