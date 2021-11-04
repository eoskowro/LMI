% Arizona State University
% MAE 598 LMIs in Control Systems

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMI for Opimal Output Feedback H_inf Control
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
X1= sdpvar(6);
Y1= sdpvar(6);
gam= sdpvar(1);
An= sdpvar(6); 
Bn= sdpvar(6,3); 
Cn= sdpvar(3,6); 
Dn= sdpvar(3,3);

% Define Matricies
mat1= [X1 eye(6);...
      eye(6) Y1];

m11= A*Y1+Y1*A'+B2*Cn+Cn'*B2'; 
m21= A'+An+(B2*Dn*C2)';
m31= (B1+B2*Dn*D21)'; 
m41= C1*Y1+D12*Cn;
m22= X1*A+A'*X1+Bn*C2+C2'*Bn';
m32= (X1*B1+Bn*D21)';
m42= C1+D12*Dn*C2;
m43= D11+D12*Dn*D21;
m33= -gam*eye(6);


mat2= [m11 m21' m31' m41';...
       m21 m22  m32' m42';...
       m31 m32  m33  m43';...
       m41 m42  m43  m33];

% Define Constraints
F = [];
F = [F, mat1>=eta*eye(12)];
F = [F, mat2<=-eta*eye(24)];

% Optimization Problem
optimize(F,gam,opt);
gam = value(gam);
disp('The predicted H-inf gain is: ');
disp(gam);

% Construct Controller
An = value(An);
Bn = value(Bn);
Cn = value(Cn);
Dn = value(Dn);
X1 = value(X1);
Y1 = value(Y1);
Y2 = eye(6); 
X2 = eye(6)-X1*Y1;

mat1 = inv([X2 X1*B2; zeros(3,6) eye(3)]);
mat2 = [An Bn; Cn Dn]- [X1*A*Y1 zeros(6,3); zeros(3,9)];
mat3 = inv([Y2' zeros(6,3); C2*Y1 eye(3)]);
mat = mat1*mat2*mat3;

% Pull out system blocks
AK2 = mat(1:6,1:6);
BK2 = mat(1:6,7:9);
CK2 = mat(7:9,1:6);
DK2 = mat(7:9,7:9);

% Form the Controller state space
DK = inv(eye(3)+DK2*D22)*DK2;
BK = BK2*(eye(3)-D22*DK);
CK = (eye(3)-DK*D22)*CK2;
AK = AK2- BK*inv(eye(3)-D22*DK)*D22*CK;
K = [AK, BK;CK, DK];
disp('The controller is:');
disp(K);

% Form CL system with Controller
Q = inv(eye(3)- D22*DK);
Acl = [A zeros(6);zeros(6) AK ] + [B2, zeros(6,3);zeros(6,3),BK]...
      *inv([eye(3), -DK; -D22, eye(3)])*[zeros(3,6), CK; C2,zeros(3,6)];
Bcl = [B1+ B2*DK*Q*D21; ...
      BK*Q*D21];
Ccl = [C1, zeros(6)]+ [D12, zeros(6,3)]*inv([eye(3), -DK;...
      -D22, eye(3)])*[zeros(3,6), CK; C2,zeros(3,6)];
Dcl = D11+D12*DK*Q*D21;

% system representation
sys = ss(Acl, Bcl, Ccl, Dcl);

% Compute Gain
hinfgain = norm(sys,inf);

% Display Results
disp('The H-inf gain with the controller is:');
disp(hinfgain);


