clc, clear all
L = 4;
h = randn(L, 1);
N = 16; 
g = [h;zeros(N-L, 1)]';
G = zeros(N, N);
epsilon = 1e-10;
for i = 1:N
    G(i,:) = circshift(g,i-1);
end
H_g = fft(g, N);
%{
F_N = fft(eye(N));
F_N = F_N/norm(F_N);
F_N_H = ifft(eye(N));
F_N_H = F_N_H/norm(F_N_H);
%}
F_N = dftmtx(N)/sqrt(N); % DFT matix
F_N_H = F_N';            % IDFT matrix
H_1 = F_N*G*F_N_H;
H_2 = F_N_H*G*F_N;

for i = 1:N
    for j = 1:N
        if abs(H_1(i, j))> epsilon
            H_1(i, j) = H_1(i, j);
        else
            H_1(i, j) = 0;
        end
    end
end

for i = 1:N
    for j = 1:N
        if abs(H_2(i, j)) > epsilon
            H_2(i, j) = H_2(i, j);
        else
            H_2(i, j) = 0;
        end
    end
end

H_1_diag = diag(H_1);
H_2_diag = diag(H_2);
diff_H_1 = H_1_diag - H_g ;
diff_H_2 = H_2_diag - H_g ;
if H_1_diag == H_g
    fprintf("true")
end
if H_1_diag == H_g
    fprintf("true")
end