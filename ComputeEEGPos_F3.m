function [F3,F3_ori,Cz,Cz_ori] = ComputeEEGPos_F3(Nz,Iz,M1,M2,msh_path)
% This function is used to get the EEG F3 locatoin using mesh2eeg.

sub_name = msh_path(end-6:end);
head_mesh = mesh_load_gmsh4([msh_path,'\',sub_name,'.msh']);
skin_mesh = mesh_extract_regions(head_mesh, 'elemtype','tri','region_idx',1005);
VerticesN = skin_mesh.nodes;
Faces = skin_mesh.triangles;

Fiducials = [Nz;Iz;M1;M2];

[~,F3,Cz] = ComputeEEGPos(Fiducials,Faces,VerticesN,2,0);

Vmin = abs(min(VerticesN(:)));
Vmax = max(VerticesN(:));

F3_ori = F3*(Vmax+Vmin)/100-Vmin;
Cz_ori = Cz*(Vmax+Vmin)/100-Vmin;

return
