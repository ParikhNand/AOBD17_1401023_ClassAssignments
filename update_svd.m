%% Function for Incremental SVD
function [U,S,V]=update_svd(U,K,V,A,B);

[m_u, n_u] = size(U);
[m_v, n_v] = size(V);
[m_a, n_a] = size(A);

m = m_u;    % Number of rows in X=U*S*V'
n = m_v;    % Number of columns in X
r = n_u;    % Rank of the given SVD
c = n_a;    % Number of columns in A and B


%%% Step 1: Find components of A in terms of left singular vectors, U.

% We decompose A into two components: M = U'*A, the component of A in U and
% Ra, the component of A orthogonal to U.  We also find P, an orthogonal basis
% of Ra. We can find these via a QR decomposition.
%
% [U,A] = Q*R = [U,P] * [I, M; 0, Ra]
% 
% The dimensions are as follows:
% U = [m_u, n_u]
% A = [m_a, n_a]

[Q,R] = qr([U,A]);

[m_q, n_q] = size(Q);
[m_r, n_r] = size(R);

% If Q and U have the same number of columns, A lies completely
% within U and there is no orthogonal component and thus we don't need
% a basis P.
%
% Otherwise, the first r = n_u columns of Q will be the same as U, and the
% remaining m_q - r columns will be the orthogonal basis P.

if (n_q == r)
    P = [];
    Ra = [];
    dimRa = 0;
else
    P = Q(:, r+1:n_q);
    Ra = R(r+1:m_r,  r+1:n_r);
    dimRa = m_r - r;
end

M = U'*A;

%%% Step 2: Find components of B in terms of right singular vectors, V.
%%% Similar to Step 1

[Q,R] = qr([V,B]);

[m_q, n_q] = size(Q);
[m_r, n_r] = size(R);

if (n_q == r)
    Q = [];
    Rb = [];
    dimRb = 0;
else
    Q = Q(:, r+1:n_q);  
    Rb = R(r+1:m_r, r+1:n_r);
    dimRb = m_r - r;
end

N = V'*B;

%%% Step 3: Construct the temporary new S and rediagonalize it

% Equation (1.4) in Brand
Kaug = [M; Ra] * [N; Rb]';
K = [K, zeros(r, dimRb); zeros(dimRa, r), zeros(dimRa, dimRb)] + Kaug;

% Now we rediagonalize via an SVD
[U2,S2,V2] = svd(K,0);

% And compute the new revised U, S, V 
U = [U,P]*U2;
S = S2;
V = [V,Q]*V2;
