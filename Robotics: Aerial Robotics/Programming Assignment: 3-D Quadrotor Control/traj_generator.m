function [ desired_state ] = traj_generator(t, state, waypoints)

persistent d d0 traj_time t_max waypoints0 w_x w_y w_z

persistent coffx coffy coffz

if nargin > 2

desired_state.vel = zeros(3,1);

desired_state.acc = zeros(3,1);

desired_state.yaw = 0;

desired_state.yawdot = 0;

d = waypoints(:,2:end) - waypoints(:,1:end-1); %Matrix 3x4, distances between points

d0 = 2 * sqrt(d(1,:).^2 + d(2,:).^2 + d(3,:).^2); %Matrix 1x4 segment times

traj_time = [0, cumsum(d0)]; %Matrix 1x5 final times at the points

waypoints0 = waypoints;

w_x = waypoints0(1,:); %Matrix 1x5 points in x direction

w_y = waypoints0(2,:); %Matrix 1x5 points in y direction

w_z = waypoints0(3,:); %Matrix 1x5 points in z direction

coffx = getCoff(w_x); %32x1

coffy = getCoff(w_y); %32x1

coffz = getCoff(w_z); %32x1

else

t_max = traj_time(end);

if(t >= t_max)

t = t_max - 0.0001;

end

t_index = find(traj_time >t,1)-1; %indeks of a segment

t_index = max(t_index,1);

scale = (t-traj_time(t_index)) / d0(t_index);

if(t == 0)

desired_state.pos = waypoints0(:,1);

desired_state.vel = zeros(3,1);

desired_state.acc = zeros(3,1);

desired_state.yaw = 0;

desired_state.yawdot = 0;

else

t0 = poly(8,0,scale)';

t1 = poly(8,1,scale)';

t2 = poly(8,2,scale)';

index = (t_index-1)*8+1:t_index*8;

desired_state.pos = [coffx(index)'*t0 ;coffy(index)'*t0 ;coffz(index)'*t0];

desired_state.vel = [coffx(index)'*t1*(1/d0(t_index)) ;coffy(index)'*t1*(1/d0(t_index)) ;coffz(index)'*t1*(1/d0(t_index))];

desired_state.acc = [coffx(index)'*t2*(1/d0(t_index)^2) ;coffy(index)'*t2*(1/d0(t_index)^2) ;coffz(index)'*t2*(1/d0(t_index)^2)];

desired_state.yaw = 0;

desired_state.yawdot = 0;

end

end %if else nargin

function [coff, A, b] = getCoff(points)

n = 4;%size(points,1)-1; % number of segments P1..n

A = zeros(8*n, 8*n);

b = zeros(1,8*n);

row = 1;

%row 1 P1(0) = W1

A(row,1:8) = polyT(8,0,0);

b(1,row) = points(1,1);

row = row+1;

%row 2 P2(0) = W2

A(row,9:16) = polyT(8,0,0);

b(1,row) = points(1,2);

row = row+1;

%row 3 P3(0) = W3

A(row,17:24) = polyT(8,0,0);

b(1,row) = points(1,3);

row = row+1;

%row 4 P4(0) = W4

A(row,25:32) = polyT(8,0,0);

b(1,row) = points(1,4);

row = row+1;

%row 5 P1(T) = W2

A(row,1:8) = polyT(8,0,1);

b(1,row) = points(1,2);

row = row+1;

%row 6 P2(T) = W3

A(row,9:16) = polyT(8,0,1);

b(1,row) = points(1,3);

row = row+1;

%row 7 P3(T) = W4

A(row,17:24) = polyT(8,0,1);

b(1,row) = points(1,4);

row = row+1;

%row8 P4(T) = W5

A(row,25:32) = polyT(8,0,1);

b(1,row) = points(1,5);

row = row+1;

%row 9 P1_dot(0) = 0

A(row, 1:8) = polyT(8,1,0);

row = row+1;

%row 10 P1_ddot(0) = 0

A(row,1:8) = polyT(8,2,0);

row = row+1;

%row 11 P1_dddot(0) = 0

A(row,1:8) = polyT(8,3,0);

row = row+1;

%row 12 P4_dot(T) = 0

A(row,25:32) = polyT(8,1,1);

row = row+1;

%row 13 P4_ddot(T) = 0

A(row,25:32) = polyT(8,2,1);

row = row+1;

%row 14 P4_dddot(T) = 0

A(row,25:32) = polyT(8,3,1);

row = row+1;

%row 15 P1_dot(T) - P2_dot(0) = 0

A(row,1:16) = [polyT(8,1,1), -polyT(8,1,0)];

row = row+1;

%row 16 P2_dot(T) - P3_dot(0) = 0

A(row,9:24) = [polyT(8,1,1), - polyT(8,1,0)];

row = row +1;

%row 17 P3_dot(T) - P4_dot(0) = 0

A(row,17:32) = [polyT(8,1,1), - polyT(8,1,0)];

row = row+1;

%row 18 P1_ddot(T) - P2_ddot(0) = 0

A(row,1:16) = [polyT(8,2,1), - polyT(8,2,0)];

row = row+1;

%row 19 P2_ddot(T) - P3_ddot(0) = 0

A(row,9:24) = [polyT(8,2,1), - polyT(8,2,0)];

row = row +1;

%row 20 P3_ddot(T) - P4_ddot(0) = 0

A(row,17:32) = [polyT(8,2,1), - polyT(8,2,0)];

row = row+1;

%row 21 P1_dot3(T) - P2_dot3(0) = 0

A(row,1:16) = [polyT(8,3,1), - polyT(8,3,0)];

row = row+1;

%row 22 P2_dot3(T) - P3_dot3(0) = 0

A(row,9:24) = [polyT(8,3,1), - polyT(8,3,0)];

row = row +1;

%row 23 P3_dot3(T) - P4_dot3(0) = 0

A(row,17:32) = [polyT(8,3,1), - polyT(8,3,0)];

row = row+1;

%row 24 P1_dot4(T) - P2_dot4(0) = 0

A(row,1:16) = [polyT(8,4,1), - polyT(8,4,0)];

row = row+1;

%row 25 P2_dot4(T) - P3_dot4(0) = 0

A(row,9:24) = [polyT(8,4,1), - polyT(8,4,0)];

row = row +1;

%row 26 P3_dot4(T) - P4_dot3(0) = 0

A(row,17:32) = [polyT(8,4,1), - polyT(8,4,0)];

row = row+1;

%row 27 P1_dot5(T) - P2_dot5(0) = 0

A(row,1:16) = [polyT(8,5,1), - polyT(8,5,0)];

row = row+1;

%row 28 P2_dot5(T) - P3_dot5(0) = 0

A(row,9:24) = [polyT(8,5,1), - polyT(8,5,0)];

row = row +1;

%row 29 P3_dot5(T) - P4_dot5(0) = 0

A(row,17:32) = [polyT(8,5,1), - polyT(8,5,0)];

row = row+1;

%row 30 P1_dot6(T) - P2_dot6(0) = 0

A(row,1:16) = [polyT(8,6,1), - polyT(8,6,0)];

row = row+1;

%row 31 P2_dot6(T) - P3_dot6(0) = 0

A(row,9:24) = [polyT(8,6,1), - polyT(8,6,0)];

row = row +1;

%row 32 P3_dot6(T) - P4_dot6(0) = 0

A(row,17:32) = [polyT(8,6,1), - polyT(8,6,0)];

%deter = det(A);

A_inv = A^(-1);

coff = A_inv*b';

function [ T ] = polyT( n, k, t)

% n is the polynom number of coefficients, k is the requested derivative and t is the actual value of t (this can be anything, not just 0 or 1).

T = zeros(n,1);

D = zeros(n,1);

% Init:

for i=1:n

D(i) = i-1;

T(i) = 1;

end

% Derivative

for j=1:k

for i=1:n

T(i) = T(i) * D(i);

if D(i) > 0

D(i) = D(i) - 1;

end

end

end

% put t value

for i=1:n

T(i) = T(i) * t^D(i);

end

T = T';

end %polyT function 1

end %getCoff function

function [ T ] = poly( n, k, t)

% n is the polynom number of coefficients, k is the requested derivative and t is the actual value of t (this can be anything, not just 0 or 1).

T = zeros(n,1);

D = zeros(n,1);

% Init:

for i=1:n

D(i) = i-1;

T(i) = 1;

end

% Derivative

for j=1:k

for i=1:n

T(i) = T(i) * D(i);

if D(i) > 0

D(i) = D(i) - 1;

end

end

end

% put t value

for i=1:n

T(i) = T(i) * t^D(i);

end

T = T';

end %polyT function 2

end %end of gen_traj