function [X,Y] = Beam_F3_pubmed(TT,NI,CF)
%TT distance from LPA to RPA
%NI distance from inion to nasion
%CF circumference
PI = 3.14159265358979323846264338327950;

R2 = TT*0.4;
R1 = NI*0.4;

A =  (R2*sind(324)+R1)/(R2*cosd(324)); %pubmed
%A = (R2*sind(324)+R1)/(R2*cosd(324));  %published
%A = (R2*sind(324)+R1/2)/R2;            %website
%A = (R2*sind(324)+R1/2)/R2*cosd(324);  %formula

B = R1*sind(288)/(R1*cosd(288)-R2);   %pubmed
%B = R1*sind(288)/(R1*cosd(288)-R2);   %published
%B = R1*sind(298)/(R1*cosd(298)-R2/2); %website
%B = R1*sind(288)/(R1*cosd(288)-R2/2); %formula

x = 1/2*(R1+R2*B)/(B-A); %pubmed
%x = 1/2*(R1+R2*B)/(B-A); %published
%x = 1/2*(R1-R2*B)/(A-B); %website
%x = 1/2*(R1-R2*B)/(A-B); %formula

y = 1/2*(A*((R1+R2*B)/(B-A))+R1); %pubmed
%y = 1/2*(A*((R1+R2*B)/(B-A))+R1); %published
%y = A*x-R1/2;                     %website
%y = A*x-R1/2;                     %formula

r = sqrt(2*x^2+2*y^2);             %pubmed
%r = sqrt(1/2*(R1+R2*B^2)/(B-A)+y^2);%published
%r = sqrt(x^2+y^2);                 %website
%r = sqrt(x^2+y^2);                 %formula

theta = atan(y/x);       %pubmed
%theta = atan((y+R1/2)/x);%published
%theta = atan(y/x);       %website
%theta = atan(y/x);       %formula

Y = r*0.9;  %pubmed
%Y = r*0.9;  %published
%Y = r*0.91; %website
%Y = r*0.9;  %formula

% X = CF*(90-theta)/360;
X = CF/4 * (90+ theta*(180/PI)) /90;

end
