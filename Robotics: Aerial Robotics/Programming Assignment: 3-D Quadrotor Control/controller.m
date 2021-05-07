function [F, M] = controller(t, state, des_state, params)
%CONTROLLER  Controller for the quadrotor
%
%   state: The current state of the robot with the following fields:
%   state.pos = [x; y; z], state.vel = [x_dot; y_dot; z_dot],
%   state.rot = [phi; theta; psi], state.omega = [p; q; r]
%
%   des_state: The desired states are:
%   des_state.pos = [x; y; z], des_state.vel = [x_dot; y_dot; z_dot],
%   des_state.acc = [x_ddot; y_ddot; z_ddot], des_state.yaw,
%   des_state.yawdot
%
%   params: robot parameters

%   Using these current and desired states, you have to compute the desired
%   controls


% =================== Your code goes here ===================

rot = state.rot;
omega = state.omega;
rd = des_state.yawdot;
rddot_T = des_state.acc;
rdot_T = des_state.vel - state.vel;
r_T = des_state.pos - state.pos;
psiT = des_state.yaw;

kp_1 = 10;
kv_1 = 20;
kp_2 = 30;
kv_2 = 500;
kp_3 = 150;
kv_3 = 500;

rddot_1d = rddot_T(1) + (kv_1*rdot_T(1)) + (kp_1*r_T(1));
rddot_2d = rddot_T(2) + (kv_2*rdot_T(2)) + (kp_2*r_T(2));
rddot_3d = rddot_T(3) + (kv_3*rdot_T(3)) + (kp_3*r_T(3));

kp_phi = 550;
kv_phi = 900;
kp_theta = 1500;
kv_theta = 1200;
kp_psi = 300;
kv_psi = 900;

phid = 1/params.gravity*((rddot_1d*sin(psiT)) - (rddot_2d*cos(psiT)));
thetad = 1/params.gravity*((rddot_1d*cos(psiT)) + (rddot_2d*sin(psiT)));
psid = psiT;

u1 = (params.mass*params.gravity) + (params.mass*rddot_3d);
u2 = [((kp_phi)*(phid-rot(1)) + (kv_phi)*(0-omega(1)));...
    ((kp_theta)*(thetad-rot(2)) + (kv_theta)*(0-omega(2)));...
    (kp_psi)*(phid-rot(3)) + (kv_psi)*(rd-omega(3))];

F = u1;
M = (params.I)*u2;



% =================== Your code ends here ===================

end
