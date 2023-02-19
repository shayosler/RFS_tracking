function [rot] = erm_b2e(eulers)
% Calculate the 3x3 rotation matrix from body frame to Earth-frame
% given Euler angles that describe the vehicle's attitude
% Inputs:
%   Eulers  Vehicle attitude euler angles [phi, theta, psi], radians
%
% Outputs:
%   Rotation matrix

phi   = eulers(1);
theta = eulers(2);
psi   = eulers(3);

rot(1,1) = cos(psi)*cos(theta);
rot(1,2) = cos(psi)*sin(theta)*sin(phi) - sin(psi)*cos(phi);
rot(1,3) = cos(psi)*sin(theta)*cos(phi) + sin(psi)*sin(phi);
rot(2,1) = sin(psi)*cos(theta);
rot(2,2) = sin(psi)*sin(theta)*sin(phi) + cos(psi)*cos(phi);
rot(2,3) = sin(psi)*sin(theta)*cos(phi) - cos(psi)*sin(phi);
rot(3,1) = -sin(theta);
rot(3,2) = cos(theta)*sin(phi);
rot(3,3) = cos(theta)*cos(phi);

end


