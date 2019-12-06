% Author: Elizabeth Skowronek
% Arizona State University
% MAE 598 LMIs in Control Systems

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bring in the data
disp('Finding the Observability grammian using an LMI');
A = rand(5,5); % Put your A matrix here
C = rand(2,5); % Put your B matrix here

% Sizes and settings
n = length(A); % Size of A matrix
m = size(C,1); % Rows of C matrix 
eps = 1*10^-8; % very small value to enforce positive definite
epsI = eps*eye(n); % matrix needed to enforce positive definite for nxn
opt = sdpsettings('solver','sedumi','verbose',0); % Optimization settings

% Define Variables
Y = sdpvar(n,n);

% Define Constraints
F = [];
F = [F Y>=epsI];
F = [F A'*Y+Y*A-C'*C<=-epsI];

% Optimization Problem
optimize(F,[],opt);

% Display Results
disp('-----------------------------------------------------------');
disp('The Observability Grammian Y is:');
disp(value(Y));

