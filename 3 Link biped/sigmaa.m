%% This function forms a q vector with give omega-
%%=sigmaa(theta1d,theta3d,dtheta1_)
function s=sigmaa(theta1d,theta3d,dtheta1_)
s=[theta1d -theta1d theta3d dtheta1_ -dtheta1_ 0];
end