function  [F3, F3_ori]= Arch_SGP_normalize_P(Nz,Iz,M1,M2,Cz,pNz,pAl,msh_path)
%====================================================================================================================
% This is function is used to get the SGP F3 based on SGP parameter.
% [input]:
%       Nz:   coordinate of Nasion.
%       Iz:   coordinate of Inion.
%       M1:   coordinate of Left Tragus.
%       M2:   coordinate of right Tragus.
%       Cz:   coordinate of Cz.
%       PNz:  distance Nasion-S point divided by distance Nasion-Inion
%       pAl:  distance M1-SGP F3 divided by distance M1-M2
%       msh_path: fullpath of .msh for one subject.
% [output] A-B: from A to B
%       F3:     SGP F3 location in nomalized space.
%       F3_ori: SGP F3 location in native space.
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
Iz_2 = roundn(Iz,-4);
Nz_2 = roundn(Nz,-4);

if EdgePts_sag_all(1,3)<EdgePts_sag_all(2,3)
    EdgePts_Nz_Iz=EdgePts_sag_all(find(sum(EdgePts_sag_all==Nz_2,2)==3):find(sum(EdgePts_sag_all==Iz_2,2)==3),:,:);
else
    EdgePts_Nz_Iz=[EdgePts_sag_all(find(sum(EdgePts_sag_all==Iz_2,2)==3):end,:);EdgePts_sag_all(1,:)];
end
EdgePts_Nz_Iz = data_smooth(Nz,Iz,M1,M2,EdgePts_Nz_Iz,1000);

NI_IZ_length = sum(sqrt(sum(diff(EdgePts_Nz_Iz).^2,2)));
rand_p = floor(size(EdgePts_Nz_Iz,1)/2);
if EdgePts_Nz_Iz(rand_p,2)>= EdgePts_Nz_Iz(rand_p+1,2)
else
    EdgePts_Nz_Iz = flipud(EdgePts_Nz_Iz);
end

cum_sum  = cumsum(sqrt(sum(diff(EdgePts_Nz_Iz).^2,2)));
inx = find(cum_sum <= pNz*NI_IZ_length);
inx = inx(end);
if cum_sum(inx)==pNz*NI_IZ_length
    S_1 = EdgePts_Nz_Iz(inx+1);
else
    dis_last = cum_sum(inx+1)- cum_sum(inx);
    last_more = pNz*NI_IZ_length - cum_sum(inx);
    y = (EdgePts_Nz_Iz(inx+2,:,:)-EdgePts_Nz_Iz(inx+1,:,:))*(last_more/dis_last);
    x = EdgePts_Nz_Iz(inx+1,:,:);
    S_1 = [x(1)+y(1),x(2)+y(2),x(3)+y(3)];
end

[EdgePts_sgp_all] = MeshPlaneIntersectPoints([M1;M2;S_1],FacesT,VerticesN);
if EdgePts_sgp_all(100,3)>EdgePts_sgp_all(1,3)
    EdgePts_sgp_all = flipud(EdgePts_sgp_all);
end

M2_2 = roundn(M2,-4);
EdgePts_M1_M2 = [EdgePts_sgp_all(find(sum(EdgePts_sgp_all==M2_2,2)==3):end,:);EdgePts_sgp_all(1,:)];
EdgePts_M1_M2 = data_smooth(Nz,Iz,M1,M2,EdgePts_M1_M2,1000);
EdgePtsS_2_Nz = sum(sqrt(sum(diff(EdgePts_M1_M2).^2,2)));
rand_p = floor(size(EdgePts_M1_M2,1)/2);
if EdgePts_M1_M2(rand_p,1)<= EdgePts_M1_M2(rand_p+1,1)
else
    EdgePts_M1_M2 = flipud(EdgePts_M1_M2);
end

cum_sum  = cumsum(sqrt(sum(diff(EdgePts_M1_M2).^2,2)));
inx = find(cum_sum <= pAl*EdgePtsS_2_Nz);
inx = inx(end);
if cum_sum(inx)==pAl*EdgePtsS_2_Nz
    S_2 = EdgePts_M1_M2(inx+1);
else
    dis_last = cum_sum(inx+1)- cum_sum(inx);
    last_more = pAl*EdgePtsS_2_Nz - cum_sum(inx);
    y = (EdgePts_M1_M2(inx+2,:,:)-EdgePts_M1_M2(inx+1,:,:))*(last_more/dis_last);
    x = EdgePts_M1_M2(inx+1,:,:);
    S_2 = [x(1)+y(1),x(2)+y(2),x(3)+y(3)];
end

F3  = S_2;
F3_ori = S_2*(Vmax+Vmin)/100-Vmin;

return
