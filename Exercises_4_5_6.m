clc
close all
clear

%% Exercise 4

syms x y e1 e2 a1 a2 b1 b2 

% USING EIGENVALUES
% x_eqn = x*(e1 -a1 * x - b1 * y);
% y_eqn = y*(e2 -a2 * y - b2 * x);
% 
% jacob_eqns = jacobian([x_eqn, y_eqn], [x ; y])
% 
% eig_val_r2 = eig(subs(jacob_eqns,[x,y],[0,e2/b2]))
% eig_val_r3 = eig(subs(jacob_eqns,[x,y],[e1/b1,0]))
% 
% [x_r4, y_r4] = solve([a1*x+b1*y==e1, a2*x+b2*y==e2], [x,y])
% eig_val_r4 = eig(subs(jacob_eqns,[x,y],[x_r4,y_r4]))
% Unsure how to determine stability with this many variables in the
% eigenvalues

% ATTEMPTED using state space matrix
A = [e1-2*a1*x-b1*y, -b1*x;-a2*y, e2-2*b2*y-a2*x]

A1 = subs(A,[x,y],[0,0]);
calc_inv_lap(A1);

A2 = subs(A,[x,y],[0,e2/b2]);
calc_inv_lap(A2)

A3 = subs(A,[x,y],[e1/b1, 0]);
calc_inv_lap(A3);

[x_r4, y_r4] = solve([a1*x+b1*y==e1, a2*x+b2*y==e2], [x,y]);
A4 = subs(A,[x,y],[x_r4,y_r4]);
calc_inv_lap(A4);

% Unsure how to determine stability with this many variables in the matrix

%% Exercise 5

syms x y a b delta gamma

x_d = a*x - b*x*y;
y_d = delta*x*y - gamma*y;

[x_coords, y_coords] = solve([x_d==0, y_d==0], [x y])

pnt1 = [x_coords(1), y_coords(1)];
pnt2 = [x_coords(2), y_coords(2)]

jacob_eqns = jacobian([x_d, y_d], [x ; y]);
eig_val_r1 = eig(subs(jacob_eqns,[x,y],[pnt1(1),pnt1(2)]));
% For fixed point (0,0), stable when alpha<0 and gamma > 0. Unstable otherwise


eig_val_r2 = eig(subs(jacob_eqns,[x,y],[pnt2(2),pnt2(2)]))
% Unsure how to determine stability with this many variables in the
% eigenvalues for fixed point (gamme/delta, alpha/beta)

%% Exercise 6

syms alpha S I R beta gamma

dS = -alpha*S*I;
dI = alpha*S*I-beta*I;
dR = beta*I;

[s_coord, i_coord, r_coord] = solve([dS==0, dI==0, dR==0], [S I R])
jacob_eqns = jacobian([dS, dI, dR], [S,I,R]);
eig_val_r1 = eig(subs(jacob_eqns,[S,I,R],[s_coord,i_coord, r_coord]));
% Fixed point (0,0,0) unstable even when R>0 due to one eigenvalue = 0

%b 
dS = -alpha*S*I + gamma*I;
dI = alpha*S*I-gamma*I;
[s_coord, i_coord] = solve([dS==0, dI==0], [S I])
fixed1 = [s_coord(1) i_coord(1)]
fixed2 = [s_coord(2) i_coord(2)]
jacob_eqns = jacobian([dS, dI], [S,I]);
eig_val_r1b = eig(subs(jacob_eqns,[S,I],fixed1))
% Fixed point (gamma/alpha,0) unstable due to two eigenvalues = 0

eig_val_r2b = eig(subs(jacob_eqns,[S,I],fixed2))
% Fixed point (0,0) unstable due to ome eigenvalues = 0


%% functions 

function inv_lap_a  = calc_inv_lap(A)
    syms s
    sI_A = [s 0;0 s]-A;
    inv_A = inv(sI_A);
    inv_lap_a = ilaplace(inv_A);
end


