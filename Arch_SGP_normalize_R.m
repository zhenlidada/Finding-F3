function [pNz,pAl] = Arch_SGP_normalize_R(Nz,Iz,M1,M2,EEGP,Cz,msh_path)
%=====================================================================================================================
% This is function is used to calculate the SGP parameters of F3 based on scalp geometry based parameter-space method.
% [input]:
%       Nz:   coordinate of Nasion.
%       Iz:   coordinate of Inion.
%       M1:   coordinate of Left Tragus.
%       M2:   coordinate of right Tragus.
%       EEGP: coordinate of EEG F3.
%       Cz:   coordinate of Cz.
%       msh_path: fullpath of .msh for one subject.
% [output] A-B: from A to B
%       PNz:  distance Nasion-S point divided bydistance Nasion-Inion
%       pAl:  distance M1-SGP F3 divided by distance M1-M2
% [Reference]
%       Jiang et al.,(2022).A scalp-measurement based parameter space: Towards locating TMS
%       coils in a clinically-friendly way.
%
% Zhen Li, 2023, Maastricht university
% zhen.li@maastrichtuniversity.nl
% zhen.li.dada@gmail.com
% =====================================================================================================================
%%
sub_name = msh_path(end-6:end);
head_mesh = mesh_load_gmsh4([msh_path,'\',sub_name,'.msh']);
skin_mesh = mesh_extract_regions(head_mesh, 'elemtype','tri','region_idx',1005);
Vertices = skin_mesh.nodes;
Faces = skin_mesh.triangles;

% Compile fiducial position matrix for input into function
Fiducials = [Nz;Iz;M1;M2];
% Normalize Vertices to range from [0-100] to facilitate computation
% and minimize error.
Vmin = abs(min(Vertices(:)));
Vmax = max(Vertices(:));
VerticesN = 100*(Vertices + Vmin)/(Vmax+Vmin);
Fiducials = 100*(Fiducials + Vmin)/(Vmax+Vmin);
Cz = 100*(Cz + Vmin)/(Vmax+Vmin);

Nz = Fiducials(1,:);
Iz = Fiducials(2,:);
M1 = Fiducials(3,:);
M2 = Fiducials(4,:);

% Convert tetrahedral mesh into triangular mesh
if size(Faces,2)==4
    TRI = TriRep(Faces,VerticesN);
    FacesT = freeBoundary(TRI);
else
    FacesT = Faces;
end

[EdgePts_sag_all] = MeshPlaneIntersectPoints([Nz;Iz;Cz],FacesT,VerticesN);
if EdgePts_sag_all(100,3)>EdgePts_sag_all(1,3)
    EdgePts_sag_all = flipud(EdgePts_sag_all);
end
Iz_2 = roundn(Iz,-4);
Nz_2 = roundn(Nz,-4);
if EdgePts_sag_all(1,3)<EdgePts_sag_all(2,3)
    EdgePts_Nz_Iz=EdgePts_sag_all(find(sum(EdgePts_sag_all==Nz_2,2)==3):find(sum(EdgePts_sag_all==Iz_2,2)==3),:,:);
else
    if sum(EdgePts_sag_all(1,:)==Nz_2,2)==3
        EdgePts_Nz_Iz=[EdgePts_sag_all(find(sum(EdgePts_sag_all==Iz_2,2)==3):end,:);EdgePts_sag_all(1,:)];
    else
        EdgePts_Nz_Iz=EdgePts_sag_all(find(sum(EdgePts_sag_all==Iz_2,2)==3):end,:);
    end
end
EdgePts_Nz_Iz = data_smooth(Nz,Iz,M1,M2,EdgePts_Nz_Iz,1000);

EEGP = 100*(EEGP + Vmin)/(Vmax+Vmin);

[EdgePts_sgp_all] = MeshPlaneIntersectPoints([M1;M2;EEGP],FacesT,VerticesN);
if EdgePts_sgp_all(100,3) > EdgePts_sgp_all(1,3)
    EdgePts_sgp_all = flipud(EdgePts_sgp_all);
end
M2_2 = roundn(M2,-4);
EdgePts_M1_M2 = [EdgePts_sgp_all(find(sum(EdgePts_sgp_all==M2_2,2)==3):end,:);EdgePts_sgp_all(1,:)];
EdgePts_M1_M2 = data_smooth(Nz,Iz,M1,M2,EdgePts_M1_M2,1000);

Fiducial = [Nz;Cz;Iz];
p1 = Fiducial(1,:);p2 = Fiducial(2,:);p3 = Fiducial(3,:);
% Compute constants for plane equation
a = (p2(1,2)-p1(1,2))*(p3(1,3)-p1(1,3))-(p3(1,2)-p1(1,2))*(p2(1,3)-p1(1,3));
b = (p2(1,3)-p1(1,3))*(p3(1,1)-p1(1,1))-(p3(1,3)-p1(1,3))*(p2(1,1)-p1(1,1));
c = (p2(1,1)-p1(1,1))*(p3(1,2)-p1(1,2))-(p3(1,1)-p1(1,1))*(p2(1,2)-p1(1,2));
d = -(a*p1(1,1)+b*p1(1,2)+c*p1(1,3)); % This would work with p1, p2, or p3
% Compute distance to all mesh vertices
Dist = (a*EdgePts_sgp_all(:,1)+b*EdgePts_sgp_all(:,2)+c*EdgePts_sgp_all(:,3)+d)/sqrt(a^2+b^2+c^2);
inx = find(Dist>0);
dot2inx1 = inx(end);
dot1inx2 = dot2inx1+1;
point2 = EdgePts_sgp_all(dot1inx2,:);
point1 = EdgePts_sgp_all(dot2inx1,:);

Fiducial = [M1;EEGP;M2];
p1 = Fiducial(1,:);p2 = Fiducial(2,:);p3 = Fiducial(3,:);
% Compute constants for plane equation
a = (p2(1,2)-p1(1,2))*(p3(1,3)-p1(1,3))-(p3(1,2)-p1(1,2))*(p2(1,3)-p1(1,3));
b = (p2(1,3)-p1(1,3))*(p3(1,1)-p1(1,1))-(p3(1,3)-p1(1,3))*(p2(1,1)-p1(1,1));
c = (p2(1,1)-p1(1,1))*(p3(1,2)-p1(1,2))-(p3(1,1)-p1(1,1))*(p2(1,2)-p1(1,2));
d = -(a*p1(1,1)+b*p1(1,2)+c*p1(1,3)); % This would work with p1, p2, or p3
% Compute distance to all mesh vertices
Dist = (a*EdgePts_sag_all(:,1)+b*EdgePts_sag_all(:,2)+c*EdgePts_sag_all(:,3)+d)/sqrt(a^2+b^2+c^2);
dot1inx3 = find(Dist<0);
dot1inx3 = dot1inx3(end);
dot2inx4 = dot1inx3+1;

point3 = EdgePts_sag_all(dot1inx3,:);
point4 = EdgePts_sag_all(dot2inx4,:);

xp1=point1(1);yp1=point1(2);zp1=point1(3);
xp2=point2(1);yp2=point2(2);zp2=point2(3);
xp3=point3(1);yp3=point3(2);zp3=point3(3);
xp4=point4(1);yp4=point4(2);zp4=point4(3);
t=((yp1-yp3)*(xp3-xp4)-(yp3-yp4)*(xp1-xp3))/((yp3-yp4)*(xp1-xp2)-(xp3-xp4)*(yp1-yp2));
CS = [xp1+(xp1-xp2)*t,yp1+(yp1-yp2)*t,zp1+(zp1-zp2)*t];
NI_IZ_length = sum(sqrt(sum(diff(EdgePts_Nz_Iz).^2,2)));
rand_p = floor(size(EdgePts_Nz_Iz,1)/2);

if EdgePts_Nz_Iz(rand_p,2)>= EdgePts_Nz_Iz(rand_p+1,2)
else
    EdgePts_Nz_Iz = flipud(EdgePts_Nz_Iz);
end
[~,inx] = min(sqrt(sum((EdgePts_Nz_Iz-point3).^2,2)));
Nz_S1_length = sum(sqrt(sum(diff(EdgePts_Nz_Iz(1:inx,:)).^2,2)))-sqrt(sum(diff([EdgePts_Nz_Iz(inx-1,:);CS]).^2,2));
pNz = Nz_S1_length/NI_IZ_length;
M1_M2_length = sum(sqrt(sum(diff(EdgePts_M1_M2).^2,2)));
rand_p = floor(size(EdgePts_M1_M2,1)/2);
if EdgePts_M1_M2(rand_p,1)<= EdgePts_M1_M2(rand_p+1,1)
else
    EdgePts_M1_M2 = flipud(EdgePts_M1_M2);
end

ins = find(EdgePts_M1_M2(:,1)<=EEGP(1));
ins = ins(end);
M1_S2more_length = sum(sqrt(sum(diff(EdgePts_M1_M2(1:ins+1,:)).^2,2)));
S2more_M1 = EdgePts_M1_M2(ins+1,:);
M1_EEGP_length = M1_S2more_length-sqrt(sum(diff([S2more_M1;EEGP]).^2,2));
pAl = M1_EEGP_length/M1_M2_length;

return