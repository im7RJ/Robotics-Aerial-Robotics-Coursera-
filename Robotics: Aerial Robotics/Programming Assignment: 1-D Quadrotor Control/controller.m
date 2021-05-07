function [ u ] = pd_controller(~, s, s_des, params)
%PD_CONTROLLER  PD controller for the height
%
%   s: 2x1 vector containing the current state [z; v_z]
%   s_des: 2x1 vector containing desired state [z; v_z]
%   params: robot parameters
u = 0;
kp = 60;
kv = 10;
e = s_des(1)-s(1);
e_dot = s_des(2)-s(2);
z_2dot_des = 0;
u = params.mass*(z_2dot_des + kp*e + kv*e_dot + params.gravity);



% FILL IN YOUR CODE HERE


end

