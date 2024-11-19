function [X,Y] = Beam_F3_published(TT,NI,CF)
%TT distance from LPA to RPA
%NI distance from inion to nasion
%CF circumference
PI = 3.14159265358979323846264338327950;

R2 = TT*0.4;
R1 = NI*0.4;

A = (R2*sind(324)+R1)/(R2*cosd(324));%published
B = R1*sind(288)/(R1*cosd(288)-R2);  %published
x = 1/2*(R1+R2*B)/(B-A);             %published
y = 1/2*(A*((R1+R2*B)/(B-A))+R1);    %published

r = sqrt(1/2*(R1+R2*B^2)/(B-A)+y^2);%published
theta = atan((y+R1/2)/x);           %published
Y = r*0.9;                          %published

X = CF/4 * (90+ theta*(180/PI)) /90;

end
