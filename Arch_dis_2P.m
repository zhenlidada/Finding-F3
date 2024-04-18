function Dis = Arch_dis_2P(P1,P2,msh_path)
%====================================================================================================================
% This function is used to the distance between to points on the head mesh.
% [input]:
%       P1:   the first point on the head mesh.
%       P2:   the second point on the head mesh.
%       msh_path: fullpath of .msh for one subject.
% [output]:
%       Dis:  the distance between to points on the head mesh.
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

Vmin = abs(min(VerticesN(:)));
Vmax = max(VerticesN(:));
VerticesN = 100*(VerticesN + Vmin)/(Vmax+Vmin);

radius = 100;
b = (P1+P2)/2;
% Get the unit vector along the direction of the line connecting (a,b)
d =(b-P2)./norm(b-P2);

% Get initial random vector which is orthogonal to the line which lies on the circle
x = b(1)/2; y = b(2)/2;
z = (-d(1)*x -d(2)*y)/d(3);

% Normalize it
p = [x,y,z];
p = p./norm(p);

% Compute a vector orthogonal to both the vector on the circle and the line connecting (a,b)
cr = cross(p,d);

% Form a system of linear equations
A = [d;p;cr];

% Just to visualize the line connecting (a,b)
% l = [linspace(0,2,100)]'*d;

% Collect points on the circle separated by 1 degree
r = zeros(360,3);

for i = 1:360
    bb = [0;cosd(i);cosd(i+90)];
    %solve for Ax =b
    x = 1*A\bb;
    x = x/norm(x);
    r(i,:) = sqrt(radius)*x; % original is x
end
node_cir = [r(:,1)+b(1),r(:,2)+b(2),r(:,3)+b(3)];

Dis_all = zeros(360,1);
P1 = 100*(P1 + Vmin)/(Vmax+Vmin);
P2 = 100*(P2 + Vmin)/(Vmax+Vmin);
node_cir = 100*(node_cir + Vmin)/(Vmax+Vmin);
for i=1:180
    [EdgePts_anytwo] = MeshPlaneIntersectPoints([P1;P2;node_cir(i,:)],Faces, VerticesN);
    dis1 = sum((EdgePts_anytwo-P1).^2,2);
    ind1 = find(dis1==min(dis1),1);
    
    dis2 = sum((EdgePts_anytwo-P2).^2,2);
    ind2 = find(dis2==min(dis2),1);
    if ind1==ind2
        P1_3 = P1*(Vmax+Vmin)/100-Vmin;
        P2_3 = P2*(Vmax+Vmin)/100-Vmin;
        cum_dis = sqrt(sum(diff([P1_3;P2_3]).^2,2));
    else
        diff1 = ind2-ind1;
        diff2 = size(EdgePts_anytwo,1)-ind2;
        if diff1< diff2
            EdgePts = EdgePts_anytwo(ind1:ind2,:);
        else
            EdgePts = [EdgePts_anytwo(ind2:end,:);EdgePts_anytwo(ind1,:)];
        end
        if sqrt(sum((P1-EdgePts(1,:)) .^ 2,2))>sqrt(sum((P2-EdgePts(1,:)) .^ 2,2))
            EdgePts = flip(EdgePts);
        end
        if sign(EdgePts(1,1)-P1(1))~=sign(EdgePts(2,1)-P1(1))||...
                sign(EdgePts(1,2)-P1(2))~=sign(EdgePts(2,2)-P1(2))||...
                sign(EdgePts(1,3)-P1(3))~=sign(EdgePts(2,3)-P1(3))
            EdgePts = [P1;EdgePts(2:end,:)];
        else
            EdgePts = [P1;EdgePts];
        end
        
        if sign(EdgePts(end-1,1)-P2(1))~=sign(EdgePts(end,1)-P2(1))||...
                sign(EdgePts(end-1,2)-P2(2))~=sign(EdgePts(end,2)-P2(2))||...
                sign(EdgePts(end-1,3)-P2(3))~=sign(EdgePts(end,3)-P2(3))
            EdgePts = [EdgePts(1:end-1,:);P2];
        else
            EdgePts = [EdgePts;P2];
        end
        EdgePts = EdgePts*(Vmax+Vmin)/100-Vmin;
        cum_dis = cumsum(sqrt(sum(diff(EdgePts).^2,2)));
    end
    Dis_all(i) = cum_dis(end);
end
ind_min = find(Dis_all==min(Dis_all(Dis_all~=0)),1);
Dis = Dis_all(ind_min);

return