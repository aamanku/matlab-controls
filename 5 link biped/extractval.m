function [alpha,tm,tp] = extractval(x)
c = [1 1 0 0.5 0];
alpha = zeros(4,7);
Ho = [zeros(4,1) eye(4)];
tm = c*[x(1) x(3) x(2) x(5) x(4)]';
tp = c*[x(1) x(2) x(3) x(4) x(5)]';
thp_dot = [x(6) x(7) x(8) x(9) x(10)]';
alpha(:,1) = [x(2) x(3) x(4) x(5)];
alpha(:,2) = Ho*[x(6); x(7); x(8); x(9); x(10)]*(tm - tp)/(c*thp_dot)/6 + alpha(:,1);
alpha(:,7) = [x(3) x(2) x(5) x(4)]';
alpha(1,3:6) = x(11:14);
alpha(2,3:6) = x(15:18);
alpha(3,3:6) = x(19:22);
alpha(4,3:6) = x(23:26);
end