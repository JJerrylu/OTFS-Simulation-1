clc, clear all
M = 128;
N = 32;
A = eye(M);
B = fft(eye(N));
A_kronecker_B = zeros(M*N,M*N);
for i = 1:M*N
    for j= 1:M*N
        A_inx_1 = ceil(i/N);
        A_inx_2 = ceil(j/N);
        B_inx_1 = mod(i, N);
        if B_inx_1 == 0
            B_inx_1 = 6;
        end
        B_inx_2 = mod(j, N);
        if B_inx_2 == 0
            B_inx_2 = 6;
        end
            A_kronecker_B(i,j)= A(A_inx_1, A_inx_2)*B(B_inx_1,B_inx_2);
    end
end