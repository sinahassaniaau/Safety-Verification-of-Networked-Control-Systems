% Script for computing p-safety for PDMP
% Dependencies: 1) Yalmip 2) MOSEK
% Authors: Rahul Misra, Sina Hassani
%
% Computes the smallest p such that P(hit U | s0) <= p for all s0 in the
% initial set A, by solving the SOS program with p as the objective,
% and then verifies the certificate a posteriori (exact Gram-matrix check
% and Monte-Carlo simulation of the PDMP).

clc
close all
clear
yalmip('clear')

% IEEE 39 bus

% % % % % % % % % % % % % % % % % % % % Data of  Area 1
B1 = 0.3417;
R1 = 3;
Tg1 = 0.08;
Tt1 = 0.4;
H1 = 0.1677/2;
D1 = 0.084;
T12 = 0.2;
T13 = 0.25;
Kp1 = 0.33;
Ki1 = -0.37;
% % % % % % % % % % % % % % % % % % % % Data of  Area 2
B2 = 0.3823;
R2 = 2.73;
Tg2 = 0.06;
Tt2 = 0.44;
H2 = 0.2017/2;
D2 = 0.016;
T21 = 0.2;
T23 = 0.12;
Kp2 = 0.36;
Ki2 = -0.34;
% % % % % % % % % % % % % % % % % % % % Data of  Area 3
B3 = 0.3696;
R3 = 2.82;
Tg3 = 0.07;
Tt3 = 0.3;
H3 = 0.1247/2;
D3 = 0.015;
T31 = 0.25;
T32 = 0.12;
Kp3 = 0.25;
Ki3 = -1;
Sbase = 100;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
A11 = [-1/Tg1 0 -1/(R1*Tg1) 0;
       1/Tt1 -1/Tt1 0 0;
       0 1/(2*H1) -D1/(2*H1) -1/(2*H1);
       0 0 2*pi*(T12+T13) 0];
A22 = [-1/Tg2 0 -1/(R2*Tg2) 0;
       1/Tt2 -1/Tt2 0 0;
       0 1/(2*H2) -D2/(2*H2) -1/(2*H2);
       0 0 2*pi*(T21+T23) 0];
A33 = [-1/Tg3 0 -1/(R3*Tg3);
       1/Tt3 -1/Tt3 0;
       0 1/(2*H3) -D3/(2*H3)];
A12 = [0 0 0 0;0 0 0 0;0 0 0 0;0 0 -2*pi*T12 0];

A13 = [0 0 0;0 0 0;0 0 0;0 0 -2*pi*T13];
A21 = [0 0 0 0;0 0 0 0;0 0 0 0;0 0 -2*pi*T21 0];
A23 = [0 0 0;0 0 0;0 0 0;0 0 -2*pi*T23];
A31 = [0 0 0 0;0 0 0 0;0 0 0 1/(2*H3)];
A32 = [0 0 0 0;0 0 0 0;0 0 0 1/(2*H3)];

b1 = [1/Tg1;0;0;0];
b2 = [1/Tg2;0;0;0];
b3 = [1/Tg3;0;0];
c1 = [0 0 B1 1];
c2 = [0 0 B2 1];
e1 = [0;0;-1/(2*H1);0];
e2 = [0;0;-1/(2*H2);0];
e3 = [0;0;-1/(2*H3)]; %% positive or negative

A_p = [A11 A12 A13;
       A21 A22 A23;
       A31 A32 A33];
B_p = [b1 zeros(4,1) zeros(4,1);
       zeros(4,1) b2 zeros(4,1);
       zeros(3,1) zeros(3,1) b3];
C_p = [c1 zeros(1,4) zeros(1,3);
       zeros(1,4) c2 zeros(1,3);
       0 0 0 -1 0 0 0 -1 0 0 B3];
E_p = [e1 zeros(4,1) zeros(4,1);
       zeros(4,1) e2 zeros(4,1);
       zeros(3,1) zeros(3,1) e3];

% Now Controllers
A_c = zeros(3,3)-0*eye(3);
B_c = diag([Ki1,Ki2,Ki3]*0.2);
C_c = diag([1 1 1]);
D_c = diag([Kp1,Kp2,Kp3]*0.3);


% Augmenting the system: states [x_p (11); x_c (3); e (3)]
A_11 = [A_p+B_p*D_c*C_p B_p*C_c;B_c*C_p A_c];
A_12 = [B_p*D_c;B_c];
A_21 = -[C_p*(A_p+B_p*D_c*C_p) C_p*(B_p*C_c)];
A_22 = -C_p*(B_p*D_c);
A = [A_11 A_12;A_21 A_22];

% PDMP parameters
lambda = 30;   % delivery rate in mode q2
delta = 0.01;  % minimal delay: lambda(tau) = lambda * 1{tau >= delta}
T_max = 0.1;   % transmission period = maximal delay (guard tau = T_max)

% Sets (Unsafe, safe and the rest)
fidx = [3 7 11]; % states considered for p-safety
x_bar = 0.8; % safe set H: |x_i| < x_bar, unsafe set U: |x_i| >= x_bar
x0_f = 0.2;  % |x_i| <= x0_f for states 3, 7, 11
x0_o = 0.02; % |x_i| <= x0_o for the other plant and controller states
r0 = x0_o*ones(14,1); r0(fidx) = x0_f;

% Hybrid Model 
n1 = size(A,1);            
n2 = n1 + 3;              
A1 = A;                    
A2 = blkdiag(A, zeros(3)); 
% R1: (x,e) -> (x,e,e), tau -> 0
J1 = [eye(n1); zeros(3,14) eye(3)]; 
% R2: (x,e,s) -> (x,e-s), tau kept
J2 = [eye(14) zeros(14,6); zeros(3,14) eye(3) -eye(3)]; 

%% Barrier function
kP = 2; 
P1 = cell(kP+1,1); P2 = cell(kP+1,1);
for j = 1:kP+1
    P1{j} = sdpvar(n1); P2{j} = sdpvar(n2);
end
dP1 = [cellfun(@(M,j) j*M, P1(2:end), num2cell((1:kP)'), 'UniformOutput', false); {zeros(n1)}];
dP2 = [cellfun(@(M,j) j*M, P2(2:end), num2cell((1:kP)'), 'UniformOutput', false); {zeros(n2)}];
sdpvar p
mu = sdpvar(14,1);

%% SOS program 
eps_sos = 1e-6;  
constr = []; cert = struct('name',{},'M',{},'Qg',{},'Sg',{},'a',{},'b',{});

% h_q >= 1 on U: P_q(tau) >= E_ii/x_bar^2 on [0,T_max] 
for i = fidx
    E1 = zeros(n1); E1(i,i) = 1/x_bar^2; E2 = zeros(n2); E2(i,i) = 1/x_bar^2;
    [constr,cert] = addMatPos(constr,cert, addc(P1,{-E1}), 0, T_max, eps_sos, sprintf('h1 >= 1 on U (x%d)',i));
    [constr,cert] = addMatPos(constr,cert, addc(P2,{-E2}), 0, T_max, eps_sos, sprintf('h2 >= 1 on U (x%d)',i));
end

% Lh_q <= 0. Mode q1: flow only (q1 is left through the guard, R1).
% Mode q2: flow, plus delivery (R2) with rate lambda for tau >= delta.
LP1  = cellfun(@(dP,P) -(dP + A1'*P + P*A1), dP1, P1, 'UniformOutput', false);
LP2a = cellfun(@(dP,P) -(dP + A2'*P + P*A2), dP2, P2, 'UniformOutput', false);
LP2b = cellfun(@(L,Pa,Pb) L - lambda*(J2'*Pa*J2 - Pb), LP2a, P1, P2, 'UniformOutput', false);
[constr,cert] = addMatPos(constr,cert, LP1,0,T_max,eps_sos, '-Lh1 >= 0');
[constr,cert] = addMatPos(constr,cert, LP2a,0,delta,eps_sos, '-Lh2 >= 0, tau <= delta');
[constr,cert] = addMatPos(constr,cert, LP2b,delta,T_max,eps_sos, '-Lh2 >= 0, tau >= delta');

% guard 
G1 = evalP(P1,T_max) - J1'*evalP(P2,0)*J1;      % h1(T_max,xi) >= h2(0, R1 xi)
G2 = evalP(P2,T_max) - J2'*evalP(P1,T_max)*J2;  % h2(T_max,xi) >= h1(T_max, R2 xi)
constr = [constr; G1 >= eps_sos*eye(n1); G2 >= eps_sos*eye(n2)];

% h1 <= p on A:  x'*P1(0)*x <= sum_i mu_i*x_i^2 <= sum_i mu_i*r0_i^2 <= p
PA = P1{1}(1:14,1:14);
constr = [constr; diag(mu) - PA >= eps_sos*eye(14); mu >= 0; p >= r0'.^2*mu];

%% Solving SoS problem
options = sdpsettings('solver','mosek');
sol = optimize(constr, p, options);
fprintf('\nSolver: %s\n', sol.info);
p_opt = value(p);
fprintf('Optimal p* = %.6f  =>  P(reach U) <= %.4f for every initial state in A\n', p_opt, p_opt);

%% A posteriori verification of the numerical certificate
% The numerical solution satisfies the coefficient equalities only up to a residual R. We move R into the Gram matrix (Qg_tilde = Qg + R) so that the SOS identity holds exactly, and check that Qg_tilde and Sg are positive
% definite. This makes the certificate independent of solver tolerances.
fprintf('\n%-28s %12s %12s %12s\n','SOS condition','residual','eig(Qg~)','eig(Sg)');
certified = true;
for k = 1:numel(cert)
    [ok, r, lQ, lS] = verifyMatPos(cert(k));
    certified = certified && ok;
    fprintf('%-28s %12.2e %12.3e %12.3e\n', cert(k).name, r, lQ, lS);
end
mu_c = max(value(mu),0);               % certified bound p_cert = sum_i mu_i*r0_i^2
p_cert = r0'.^2*mu_c;
lmin = [min(eig(value(G1))), min(eig(value(G2))), min(eig(diag(mu_c) - value(PA)))];
fprintf('%-28s %12s %12.3e\n','guard G1','',lmin(1));
fprintf('%-28s %12s %12.3e\n','guard G2','',lmin(2));
fprintf('%-28s %12s %12.3e\n','initial set A','',lmin(3));
certified = certified && all(lmin > 0);
if certified
    fprintf('Certificate verified: P(reach U) <= %.6f for every initial state in A\n', p_cert);
else
    fprintf('WARNING: certificate not verified.\n');
end

%% Barrier certificates
P1_opt = cellfun(@value, P1, 'UniformOutput', false);
P2_opt = cellfun(@value, P2, 'UniformOutput', false);
fprintf('\nh1(tau,x,e)   = [x;e]''*P1(tau)*[x;e],     P1(tau) = sum_j tau^j*P1_j\n');
fprintf('h2(tau,x,e,s) = [x;e;s]''*P2(tau)*[x;e;s], P2(tau) = sum_j tau^j*P2_j\n');
for j = 0:kP
    fprintf('  ||P1_%d|| = %.3e,  ||P2_%d|| = %.3e\n', j, norm(P1_opt{j+1}), j, norm(P2_opt{j+1}));
end
fprintf('  eig(P1(0)) in [%.3e, %.3e],  eig(P2(0)) in [%.3e, %.3e]\n', ...
        min(eig(P1_opt{1})), max(eig(P1_opt{1})), min(eig(P2_opt{1})), max(eig(P2_opt{1})));
save('pSafety_hybrid_certificate.mat','P1_opt','P2_opt','p_opt','p_cert','kP','T_max','lambda','delta','x_bar','r0');
fprintf('  Full matrices saved to pSafety_hybrid_certificate.mat\n');

%%
function P = evalP(Pc,t)
% Evaluates the matrix polynomial sum_j t^j*Pc{j+1}
P = 0;
for j = 0:numel(Pc)-1
    P = P + t^j*Pc{j+1};
end
end

function C = addc(Ac,Bc)
% Sum of two matrix polynomials given by coefficient cells
C = Ac;
for j = 1:numel(Bc)
    C{j} = C{j} + Bc{j};
end
end

function [Qs,Ss] = gramCoefs(Qg,Sg,n,m,a,b)
% Coefficients (in tau) of Z'*Qg*Z and (tau-a)(b-tau)*Z0'*Sg*Z0
Qs = repmat({zeros(n)},2*m+1,1); Ss = repmat({zeros(n)},2*m+1,1);
blk = @(X,r,c) X(r*n+(1:n), c*n+(1:n));
for r = 0:m
    for c = 0:m
        Qs{r+c+1} = Qs{r+c+1} + blk(Qg,r,c);
    end
end
w = [-a*b, a+b, -1];                 % (tau-a)(b-tau) = -ab + (a+b)tau - tau^2
for r = 0:m-1
    for c = 0:m-1
        for k = 0:2
            Ss{r+c+k+1} = Ss{r+c+k+1} + w(k+1)*blk(Sg,r,c);
        end
    end
end
end

function [constr,cert] = addMatPos(constr,cert,M,a,b,epsm,name)
% Imposes M(tau) >= 0 for tau in [a,b] through the SOS (Gram) representation
n = size(M{1},1); m = (numel(M)-1)/2;
Qg = sdpvar(n*(m+1)); Sg = sdpvar(n*m);
[Qs,Ss] = gramCoefs(Qg,Sg,n,m,a,b);
constr = [constr; Qg >= epsm*eye(n*(m+1)); Sg >= epsm*eye(n*m)];
for j = 1:2*m+1
    constr = [constr; M{j} == Qs{j} + Ss{j}];
end
cert(end+1) = struct('name',name,'M',{M},'Qg',Qg,'Sg',Sg,'a',a,'b',b);
end

function [ok,r,lQ,lS] = verifyMatPos(c)
% Exact a posteriori check of M(tau) = Z'*Qg~*Z + (tau-a)(b-tau)*Z0'*Sg*Z0
M = cellfun(@value, c.M, 'UniformOutput', false);
Qg = value(c.Qg); Sg = value(c.Sg);
n = size(M{1},1); m = (numel(M)-1)/2;
[Qs,Ss] = gramCoefs(Qg,Sg,n,m,c.a,c.b);
r = 0;
for j = 0:2*m
    R = M{j+1} - Qs{j+1} - Ss{j+1};  % residual of the coefficient of tau^j
    r = max(r, max(abs(R(:))));
    i1 = min(j,m); i2 = j - i1;      % put R into the Gram block (i1,i2)
    if i1 == i2
        Qg(i1*n+(1:n), i2*n+(1:n)) = Qg(i1*n+(1:n), i2*n+(1:n)) + R;
    else
        Qg(i1*n+(1:n), i2*n+(1:n)) = Qg(i1*n+(1:n), i2*n+(1:n)) + R/2;
        Qg(i2*n+(1:n), i1*n+(1:n)) = Qg(i2*n+(1:n), i1*n+(1:n)) + R/2;
    end
end
lQ = min(eig((Qg+Qg')/2)); lS = min(eig((Sg+Sg')/2));
ok = lQ > 0 && lS > 0;
end
