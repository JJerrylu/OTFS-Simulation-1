clc, clear all
%% part 1
P1 = [0 1 0 0 0 0].';

P  = zeros(length(P1), length(P1));
for i = 1: length(P1)
P(:, i) = circshift(P1, i-1);
end
P_2 = P*P;
%% part 2
x = 1:15;
x = x.'; 
M = 5;
N = 3;
MN = M*N;

E = zeros(M, N, MN);
P_matrix = zeros(MN, MN);
% P_matrix = [1 0 0 0 0 0; 0 0 1 0 0 0 ; 0 0 0 0 1 0; ...
%            0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1]; 
% x_permuted = P_matrix* x;
% x_permuted_2 = P_matrix.'*x_permuted;
inx = 1;
j = 1;
for i = 1:MN
    inx_1 = inx + ceil(i/M)-1;
    inx_2 = j +inx_1 -1; 
    P_matrix(i, inx_2) = 1;
    j = j + N;
    j = mod(j, MN);
end
x_permuted = P_matrix* x;
x_permuted_2 = P_matrix.'*x_permuted;


