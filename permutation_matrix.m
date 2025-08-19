clc, clear all
P1 = [0 1 0 0 0 0].';

P  = zeros(length(P1), length(P1));
for i = 1: length(P1)
P(:, i) = circshift(P1, i-1);
end
P_2 = P*P;