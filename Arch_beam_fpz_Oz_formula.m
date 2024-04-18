function [F3,F3_ori,X,Y,CM_ori,NI_ori,TT_ori,Cz] = Arch_beam_fpz_Oz_formula(Nz,Iz,M1,M2,Cz,msh_path)
%====================================================================================================================
% This is function is used to get the SGP F3 based on the SGP parameters.
% [input]:
%       Nz:   coordinate of Nasion.
%       Iz:   coordinate of Inion.
%       M1:   coordinate of Left Tragus.
%       M2:   coordinate of right Tragus.
%       Cz:   coordinate of Cz.
%       PNz:  ratio of distance Nasion-S point divided by distance Nasion-Inion
%       pAl:  ratio of distance M1-SGP F3 divided by distance M1-M2
%       msh_path: fullpath of .msh for one subject.
% [output]:
%       F3:     SGP F3 location in the nomalized space.
%       F3_ori: SGP F3 location in the native space.
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
VerticesN = skin_mesh.nodes;
Faces = skin_mesh.triangles;

Fiducials = [Nz;Iz;M1;M2];

Vmin = abs(min(VerticesN(:)));
Vmax = max(VerticesN(:));
VerticesN = 100*(VerticesN + Vmin)/(Vmax+Vmin);
Fiducials = 100*(Fiducials + Vmin)/(Vmax+Vmin);
Cz = 100*(Cz + Vmin)/(Vmax+Vmin);

if size(Faces,2)==4
    TRI = TriRep(Faces,VerticesN);
    FacesT = freeBoundary(TRI);
else
    FacesT = Faces;
end

Nz = Fiducials(1,:);
Iz = Fiducials(2,:);
M1 = Fiducials(3,:);
M2 = Fiducials(4,:);

[EdgePts_sag_all] = MeshPlaneIntersectPoints([Nz;Iz;Cz], Faces, VerticesN);
plot3(EdgePts_sag_all(:,1),EdgePts_sag_all(:,2),EdgePts_sag_all(:,3),'.','color','b')
Iz_2 = roundn(Iz,-4);
Nz_2 = roundn(Nz,-4);
if EdgePts_sag_all(1,3)<EdgePts_sag_all(2,3)
    EdgePts = EdgePts_sag_all(find(sum(EdgePts_sag_all==Nz_2,2)==3):find(sum(EdgePts_sag_all==Iz_2,2)==3),:,:);
    EdgePts = flipud(EdgePts);
else
    EdgePts=[EdgePts_sag_all(find(sum(EdgePts_sag_all==Iz_2,2)==3):end,:);EdgePts_sag_all(1,:)];
end
EdgePts = data_smooth(Nz,Iz,M1,M2,EdgePts,1000);
plot3(EdgePts(:,1),EdgePts(:,2),EdgePts(:,3),'.','color','r')
NI_ori = cumsum(sqrt(sum(diff(EdgePts*(Vmax+Vmin)/100-Vmin).^2,2)));
NI_ori = NI_ori(end);

EdgePts = EdgePts*(Vmax+Vmin)/100-Vmin;

X = sum(sqrt(sum(diff(EdgePts).^2,2)));
cum_sum  = cumsum(sqrt(sum(diff(EdgePts).^2,2)));
inx = find(cum_sum <= 0.1*X);
inx = inx(end);
if cum_sum(inx)==0.1*X
    headpoint = EdgePts(inx+1);
else
    dis_last = cum_sum(inx+1)- cum_sum(inx);
    last_more = 0.1*X - cum_sum(inx);
    y = (EdgePts(inx+2,:,:)-EdgePts(inx+1,:,:))*(last_more/dis_last);
    x = EdgePts(inx+1,:,:);
    headpoint = [x(1)+y(1),x(2)+y(2),x(3)+y(3)];
end

if headpoint(1,2)< 0
    backhead_end = headpoint;
else
    forehead_end = headpoint;
end

EdgePts = flipud(EdgePts);
X = sum(sqrt(sum(diff(EdgePts).^2,2)));
cum_sum  = cumsum(sqrt(sum(diff(EdgePts).^2,2)));
inx = find(cum_sum <= 0.1*X);
inx = inx(end);
if cum_sum(inx)==(0.1*X)
    headpoint = EdgePts(inx+1);
else
    dis_last = cum_sum(inx+1)- cum_sum(inx);
    last_more = 0.1*X - cum_sum(inx);
    y = (EdgePts(inx+2,:,:)-EdgePts(inx+1,:,:))*(last_more/dis_last);
    x = EdgePts(inx+1,:,:);
    headpoint = [x(1)+y(1),x(2)+y(2),x(3)+y(3)];
end

if headpoint(1,2)<0
    backhead_end = headpoint;
else
    forehead_end = headpoint;
end

forehead_end = 100*(forehead_end + Vmin)/(Vmax+Vmin);
backhead_end = 100*(backhead_end + Vmin)/(Vmax+Vmin);

p1 = forehead_end;
p2 = backhead_end;
p3 = Cz;

% Compute constants for plane equation
a = (p2(1,2)-p1(1,2))*(p3(1,3)-p1(1,3))-(p3(1,2)-p1(1,2))*(p2(1,3)-p1(1,3));
b = (p2(1,3)-p1(1,3))*(p3(1,1)-p1(1,1))-(p3(1,3)-p1(1,3))*(p2(1,1)-p1(1,1));
c = (p2(1,1)-p1(1,1))*(p3(1,2)-p1(1,2))-(p3(1,1)-p1(1,1))*(p2(1,2)-p1(1,2));
d = -(a*p1(1,1)+b*p1(1,2)+c*p1(1,3)); % This would work with p1, p2, or p3

x=p2(1,1)-p1(1,1);
y=p2(1,2)-p1(1,2);
z=p2(1,3)-p1(1,3);
a2=b*z-c*y ;
b2=c*x-a*z;
c2=a*y-b*x;
d2 = -(a2*p1(1,1)+b2*p1(1,2)+c2*p1(1,3));

[EdgePts_coro_all_head] = MeshPlaneIntersectPoints2(a2,b2,c2,d2,[forehead_end;backhead_end],Faces, VerticesN);
EdgePts_coro_all_head = data_smooth(Nz,Iz,M1,M2,EdgePts_coro_all_head,1000);
CM_ori = cumsum(sqrt(sum(diff(EdgePts_coro_all_head*(Vmax+Vmin)/100-Vmin).^2,2)));
CM_ori = CM_ori(end);

[EdgePts_M1M2] = MeshPlaneIntersectPoints([M1;M2;Cz],Faces, VerticesN);
M1_2 = roundn(M1,-4);
M2_2 = roundn(M2,-4);
if EdgePts_M1M2(1,3)<EdgePts_M1M2(2,3)
    EdgePts_M1M2 = flipud(EdgePts_M1M2);
end
if sum(EdgePts_M1M2(1,:)==M1_2,2)==3
    EdgePts = [EdgePts_M1M2(find(sum(EdgePts_M1M2==M2_2,2)==3):end,:);EdgePts_M1M2(1,:)];
else
    EdgePts = EdgePts_M1M2(find(sum(EdgePts_M1M2==M2_2,2)==3):end,:);
end
plot3(EdgePts_M1M2(:,1),EdgePts_M1M2(:,2),EdgePts_M1M2(:,3),'.','color','b')

[EdgePts] = Fix_TT2(EdgePts);
plot3(EdgePts(:,1),EdgePts(:,2),EdgePts(:,3),'.','color','g')
EdgePts = data_smooth(Nz,Iz,M1,M2,EdgePts,1000);
plot3(EdgePts(:,1),EdgePts(:,2),EdgePts(:,3),'.','color','r')
TT_ori = cumsum(sqrt(sum(diff(EdgePts*(Vmax+Vmin)/100-Vmin).^2,2)));
TT_ori = TT_ori(end);

[X,Y] = Beam_F3_formula(TT_ori,NI_ori,CM_ori);

if EdgePts_coro_all_head(1,1)< EdgePts_coro_all_head(2,1)
    EdgePts_coro_all_head = flipud(EdgePts_coro_all_head);
end
if forehead_end(1,1)<=EdgePts_coro_all_head(1,1)
    EdgePts = [forehead_end;EdgePts_coro_all_head(2:end,:)];
else
    EdgePts = [forehead_end;EdgePts_coro_all_head];
end

EdgePts = EdgePts*(Vmax+Vmin)/100-Vmin;
cum_sum  = cumsum(sqrt(sum(diff(EdgePts).^2,2)));
inx = find(cum_sum <= X);
inx = inx(end);
if cum_sum(inx)==X
    Xpoint = EdgePts(inx+1);
else
    dis_last = cum_sum(inx+1)- cum_sum(inx);
    last_more = X - cum_sum(inx);
    y = (EdgePts(inx+2,:)-EdgePts(inx+1,:))*(last_more/dis_last);
    x = EdgePts(inx+1,:);
    Xpoint = [x(1)+y(1),x(2)+y(2),x(3)+y(3)];
end

if EdgePts(1,1)< EdgePts(2,1)
    EdgePts = flipud(EdgePts);
end

Xpoint = 100*(Xpoint + Vmin)/(Vmax+Vmin);
[EdgePts_sf_all_head] = MeshPlaneIntersectPoints([[Cz(1),Cz(2),Cz(3)+20];[Cz;Xpoint]],Faces, VerticesN);

if EdgePts_sf_all_head(1,2)<EdgePts_sf_all_head(2,2)
else
    EdgePts_sf_all_head = flipud(EdgePts_sf_all_head);
    EdgePts_sf_all_head = [EdgePts_sf_all_head(end,:);EdgePts_sf_all_head(1:end-1,:)];
end
indcz = find(EdgePts_sf_all_head(:,2)>=Cz(1,2)&EdgePts_sf_all_head(:,1)<=Cz(1,1));
EdgePts_sf_all_head = EdgePts_sf_all_head(indcz,:);

if EdgePts_sf_all_head(indcz(1),2)<=Cz(1,2)
    EdgePts_sf_all_head = [Cz;EdgePts_sf_all_head(2:end,:)];
else
    EdgePts_sf_all_head = [Cz;EdgePts_sf_all_head];
end

EdgePts_sf_all_head = EdgePts_sf_all_head*(Vmax+Vmin)/100-Vmin;
cum_sum  = cumsum(sqrt(sum(diff(EdgePts_sf_all_head).^2,2)));
inx = find(cum_sum <= Y);
inx = inx(end);
if cum_sum(inx)==Y
    inx_Nz = inx+1;
    Ypoint = EdgePts_sf_all_head(inx_Nz);
    Ypoint_all = Ypoint;
else
    dis_last = cum_sum(inx+1)- cum_sum(inx);
    last_more = Y - cum_sum(inx);
    y = (EdgePts_sf_all_head(inx+2,:,:)-EdgePts_sf_all_head(inx+1,:,:))*(last_more/dis_last);
    x = EdgePts_sf_all_head(inx+1,:,:);
    Ypoint = [x(1)+y(1),x(2)+y(2),x(3)+y(3)];
    Ypoint_all = Ypoint;
end

F3_ori = Ypoint;
F3 = 100*(Ypoint_all + Vmin)/(Vmax+Vmin);
Cz = Cz*(Vmax+Vmin)/100-Vmin;

return

