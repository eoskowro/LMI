% Author: Elizabeth Skowronek
% Arizona State University
% MAE 598 LMIs in Control Systems

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bring in the data
disp('Finding the Controllability Grammian using an LMI');
A = rand(5,5); % Put your A matrix here
B = rand(5,2); % Put your B matrix here

% Sizes and settings
n = length(A); % Size of A matrix
m = size(B,2); % Columns of B matrix 
eps = 1*10^-8; % very small value to enforce positive definite
epsI = eps*eye(n); % matrix needed to enforce positive definite for nxn
opt = sdpsettings('solver','sedumi','verbose',0); % Optimization settings

% Define variables
W = sdpvar(n,n);

% Define Constraints
F = [];
F = [F W>=epsI]; 
F = [F A*W+W*A'-B*B'<=-epsI];

% Optimization Problem
optimize(F,[],opt); 

% Display Results
disp('The controllability grammian W is:');
disp(value(W));

