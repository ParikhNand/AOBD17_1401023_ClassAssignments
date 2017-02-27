% AOBD class - Online SVD
% Created & Written by Nand Parikh

clc; clear all;

disp('Welcome to Online SVD Program!\nFirst create m*n matrix!\n')
m = input('Enter dim. m:');
n = input('Enter dim. n:');

X = randn(m,n);         % Data Matrix
[U,S,V] = svd(X,0);     % SVD Decomposition
r = randn(n,1);         % Append a row to X
U = [U; zeros(1,n)];    
X0 = [X; zeros(1,n)];

% We want to find SVD of [X; r'] = X0 + A*B'
A = zeros(m+1,1); A(m+1,1) = 1;
B = r;

[U2,S2,V2] = update_svd(U,S,V,A,B);
disp('Success!');

new_X = X0 + A*B';          % To cross-check
[AU,AS,AV] = svd(new_X,0);

disp(S2);       % Using MatLAB SVD
disp(AS);       % Using function