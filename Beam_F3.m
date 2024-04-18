function [X,Y] = Beam_F3(TT,NI,CF)
%====================================================================================================================
% This is function is used to calculate X and Y based on Beam F3 method - webpage version.
% [input]:
%       TT, distance from Left tragus to Right tragus
%       NI, distance from nasion to Inion
%       CF, circumference
% [output] A-B: from A to B
%       X:  the distance X from FPz to point A at the FPz-Oz level along the head circumference.
%       Y:  the distance Y from the Cz position to a specific point, represented derived F3 position.
% [Reference]
%       William Beam et al.,(2009).An efficient and accurate new method for locating the F3 position
% for prefrontal TMS applications.
% ====================================================================================================================
% Zhen Li, 2023, Maastricht university
% zhen.li@maastrichtuniversity.nl
% zhen.li.dada@gmail.com
% =====================================================================================================================
PI = 3.14159265358979323846264338327950;

R2 = TT*0.4;
R1 = NI*0.4;

A = (R2*sind(324)+R1/2)/R2;
B = R1*sind(298)/(R1*cosd(298)-R2/2);
x = 1/2*(R1-R2*B)/(A-B);
y = A*x-R1/2;
r = sqrt(x^2+y^2);
theta = atan(y/x);
Y = r*.91;

X = CF/4 * (90+ theta*(180/PI)) /90;

end
