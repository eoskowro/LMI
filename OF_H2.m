% Arizona State University
% MAE 598 LMIs in Control Systems

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMI for Opimal Output Feedback H_2 Control
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
X1 = sdpvar(6);
Y1 = sdpvar(6);
Z = sdpvar(6);
gam = sdpvar(1);
An = sdpvar(6,6,'full'); 
Bn = sdpvar(6,3,'full');
Cn = sdpvar(3,6,'full'); 
Dn = sdpvar(3,3,'full');

% Define Matricies
m11 = A*Y1 + Y1*A' + B2*Cn + (Cn')*B2'; 
m21 = A'+An+(B2*Dn*C2)';
m31 = (B1+B2*Dn*D21)';
m22 = X1*A + (A')*X1 + Bn*C2 + (C2')*Bn';
m32 = (X1*B1+Bn*D21)';
m33 = -eye(6);
mm11 = Y1;
mm21 = eye(6);
mm22 = X1;
mm31 = C1*Y1+D12*Cn;
mm32 = C1+D12*Dn*C2;
mm33 = Z;

mat1=[m11 m21' m31';...
      m21 m22  m32';...
      m31 m32  m33];
mat2 = [mm11 mm21' mm31';...
      mm21 mm22  mm32';...
      mm31 mm32  mm33];  
mat3 = trace(Z);
mat4 = D11+D12*Dn*D21;

% Define Constraints
F = [];
F = [F, mat1<=eta*eye(18)];
F = [F, mat2>=eta*eye(18)];
F = [F, mat3<=gam];
F = [F, mat4==0];

% Optimization Problem
sol= optimize(F,gam,opt);
gam = value(sqrt(gam));

% Display Results
disp('The predicted H-2 gain is: ');
disp(gam);

% Bring in Values
An = value(An);
Bn = value(Bn);
Cn = value(Cn);
Dn = value(Dn);
X1 = value(X1);
Y1 = value(Y1);
Y2 = eye(6); 
X2 = eye(6)-X1*Y1;

% Define Matricies
mat1 = inv([X2 X1*B2; zeros(3,6) eye(3)]);
mat2 = ([An Bn; Cn Dn]- [X1*A*Y1 zeros(6,3); zeros(3,9)]);
mat3 = inv([Y2' zeros(6,3); C2*Y1 eye(3)]);
mat = mat1*mat2*mat3;

% Define the K2 system
AK2 = mat(1:6,1:6);
BK2 = mat(1:6,7:9);
CK2 = mat(7:9,1:6);
DK2 = mat(7:9,7:9);

% Form the controller
DK = inv(eye(3)+DK2*D22)*DK2;
BK = BK2*(eye(3)-D22*DK);
CK = (eye(3)-DK*D22)*CK2;
AK = AK2- BK*inv(eye(3)-D22*DK)*D22*CK;

K = [AK, BK;CK, DK];
disp('The controller is:');
disp(K);

% Form the CL system
Q = inv(eye(3)- D22*DK);
Acl = [A zeros(6);zeros(6) AK ] + [B2, zeros(6,3);zeros(6,3),BK]*...
       inv([eye(3), -DK; -D22, eye(3)])*[zeros(3,6), CK; C2,zeros(3,6)];
Bcl = [B1+ B2*DK*Q*D21;...
       BK*Q*D21];
Ccl = [C1, zeros(6)]+ [D12, zeros(6,3)]*inv([eye(3), -DK; -D22, eye(3)])...
       *[zeros(3,6), CK; C2,zeros(3,6)];

% Define the system state space
sys = ss(Acl,Bcl,Ccl,0);

% Compute the norm
H2norm = norm(sys,2);

% Display Results
disp('The H-2 gain with the controller is:');
disp(H2norm)