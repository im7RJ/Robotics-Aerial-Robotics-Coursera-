function [ u1, u2 ] = controller(~, state, des_state, params)
%CONTROLLER  Controller for the planar quadrotor
%
%   state: The current state of the robot with the following fields:
%   state.pos = [y; z], state.vel = [y_dot; z_dot], state.rot = [phi],
%   state.omega = [phi_dot]
%
%   des_state: The desired states are:
%   des_state.pos = [y; z], des_state.vel = [y_dot; z_dot], des_state.acc =
%   [y_ddot; z_ddot]
%
%   params: robot parameters

%   Using these current and desired states, you have to compute the desired
%   controls

z_T_dot_dot = des_state.acc(2);
phi_c_dot_dot = 0;
phi_c_dot = 0;
y_T_dot_dot = des_state.acc(1);
y_T_dot = des_state.vel(1);
z_T_dot = des_state.vel(2);
y_T = des_state.pos(1);
z_T = des_state.pos(2);
y_dot = state.vel(1);
z_dot = state.vel(2);
y = state.pos(1);
z = state.pos(2);
phi_dot = state.omega;
phi = state.rot;

kp_y = 5;
kv_y = 10;
kp_z = 300;
kv_z = 130;
kp_phi = 800;
kv_phi = 100;

ep_y = des_state.pos(1) - state.pos(1);
ev_y = des_state.vel(1) - state.vel(1);
ep_z = des_state.pos(2) - state.pos(2);
ev_z = des_state.vel(2) - state.vel(2);
edot_phi = phi_c_dot - state.omega;

phi_c = (-1/params.gravity)*(y_T_dot_dot+kv_y*(y_T_dot-y_dot)+kp_y*(y_T-y));
u1 = params.mass*(params.gravity+(z_T_dot_dot+kv_z*(z_T_dot-z_dot)+kp_z*(z_T-z)));
u2 = params.Ixx*(kv_phi*(-phi_dot)+kp_phi*(phi_c-phi));

% FILL IN YOUR CODE HERE

end

