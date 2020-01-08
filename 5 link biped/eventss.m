function [value,isterminal,direction] = eventss(t,q)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
global g lT lf lt MT Mf Mt IT If It pMT pMf pMt thetap thetam a
value=lt*cos(pi-q(1)-q(2)-q(4))+lf*cos(pi-q(1)-q(2))-lf*cos(q(1)+q(3)-pi)-lt*cos(q(1)+q(3)+q(5)-pi);
isterminal=1;
direction=0;
end