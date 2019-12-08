
% Author: Elizabeth Skowronek
% Arizona State University
% MAE 598 LMIs in Control Systems

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMI for Opimal FSF H_inf Control
% Bring in Data
A = [-1 1 0 1 0 1;...
    -1 -2 -1 0 0 1;...
    1 0 -2 -1 1 1;...
    -1 1 -1 -2 0 0;...
    -1 -1 1 1 -2 -1;...
    0 -1 0 0 -1 -3];
B = [0 -1 -1;...
    0 0 0;...
    -1 1 1;...
    -1 0 0;...
    0 0 1;...
    -1 1 1];
C = [0 1 0 -1 -1 -1;...
    0 0 0 -1 0 0;...
    1 0 0 0 -1 0];
D = zeros(3);

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
Y = sdpvar(6);
Z = sdpvar(3,6);
gam = sdpvar(1);

% Define matricies
mat = [Y*A'+A*Y+Z'*B2'+B2*Z B1 Y*C1'+Z'*D12';...
      B1' -gam*eye(6) D11';...
      C1*Y+D12*Z D11 -gam*eye(6)];
  
% Define Constraints
F = [];
F = [F, Y>=eta*eye(6)];
F = [F, mat<=-eta*eye(18)];

% Optimization Problem
optimize(F,gam,opt);
Y = value(Y);
Z = value(Z);
F = Z*inv(Y); % controller
gam = value(gam);

% Display Results
disp('The controller F is: ');
disp(F);
disp('The H-inf gain is:');
disp(gam);