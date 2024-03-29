% Arizona State University
% MAE 598 LMIs in Control Systems

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMI for Opimal FSF H_2 Control
% Bring in Data
A = [];
B = [];
C = [];
D = [];

% Determine sizes
ns = size(A,1);   % number of states
nai = size(B,2);  % number of actuator inputs
nmo = size(C,1);  % number of measured outputs
nd = nai+nmo;     % number of disturbances
nro = nmo+nai;    % number of regulated outputs
eta = 0.0001;

% 9 Matrix Representation
B1 = [B zeros(ns,nmo)];
B2 = B;
C1 = [C;...
     zeros(nai,ns)];
C2 = C;
D11 = [D zeros(nai);...
      zeros(nmo) zeros(nmo)];
D12 = [D;...
      eye(nmo)];
D21 = [D eye(nmo)];
D22 = D;

% Settings
opt = sdpsettings('verbose',0,'solver','sedumi');

% Define Variables
gam= sdpvar(1);
X = sdpvar(6);
W = sdpvar(6);
Z = sdpvar(3,6);

% Define matricies
mat1 = [A B2]*[X;Z] + [X Z']*[A';B2']+ B1*B1';
mat2 = [X (C1*X+D12*Z)';C1*X+D12*Z W];
mat3 = trace(W);

% Define Constraints
F = [];
F = [F, X>=eta*eye(6)];
F = [F, mat1<=-eta*eye(6)];
F = [F, mat2>=eta*eye(12)];
F = [F, mat3<=gam];

% Optimization Problem
optimize(F,gam,opt);
gam= sqrt(gam);
disp('The H-2 gain is: ');
disp(value(gam));
K= value(Z)*inv(value(X));

