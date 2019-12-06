% Author: Elizabeth Skowronek
% Arizona State University
% MAE 598 LMIs in Control Systems

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bring in the data
disp('Implementation of the dual KYP lemma to find an H-inf bound');
A = rand(5,5); % Put your A matrix here
B = rand(5,2); % Put your B matrix here
C = rand(2,5); % Put your C matrix here
D = zeros(2,2); % Put your D matrix here

% Sizes and settings
n = length(A); % Size of A matrix
m = size(B,2); % Columns of B matrix 
eps = 1*10^-8; % very small value to enforce positive definite
epsI = eps*eye(n); % matrix needed to enforce positive definite for nxn
opt = sdpsettings('solver','sedumi','verbose',0); % Optimization settings

% Define variables
X = sdpvar(n,n);
gam = sdpvar(1);

% Define matricies
mat = [A'*X+X*A, X*B, C'; B'*X, -gam*eye(m), D'; C, D, -gam*eye(m)];

% Define Constraints
F = [];
F = [F X>=epsI]; 
F = [F mat<=-eps*eye(size(mat))];

% Optimization Problem
optimize(F,gam,opt); 

% Display Results
disp('The value of gamma is:');
disp(value(gam));

