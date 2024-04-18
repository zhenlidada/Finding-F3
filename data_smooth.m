function EdgePts_P = data_smooth(Nz,Iz,M1,M2,EdgePts,iteration)
%====================================================================================================================
% This function is used to smooth the arch to avoid overfitting, and to be more realistic.
% [input]:
%       Nz:   coordinate of Nasion.
%       Iz:   coordinate of Inion.
%       M1:   coordinate of Left Tragus.
%       M2:   coordinate of right Tragus.
%       EdgePts:    a series of points created by the intersection plane of three points and a head mesh
%       iterations: iterations times during data smoothing.
% [output]:
%       EdgePts_P: smoothed
%
% Zhen Li, 2023, Maastricht university
% zhen.li@maastrichtuniversity.nl
% zhen.li.dada@gmail.com
% =====================================================================================================================
%%
midp = (Nz+Iz+M1+M2)/4;
midp = [midp(1),midp(2),midp(3)+30];
dis_cz_ori = sqrt(sum((midp-EdgePts) .^ 2,2));

EdgePts_P1= EdgePts(1,:);
EdgePts_P2 = (EdgePts(1:end-2,:)+EdgePts(2:end-1,:)+EdgePts(3:end,:))/3;
EdgePts_P3 = EdgePts(end,:);
EdgePts_P = [EdgePts_P1;EdgePts_P2;EdgePts_P3];

dis_midp_ir = sqrt(sum((midp-EdgePts_P) .^ 2,2));
indreplace = find(dis_midp_ir>=dis_cz_ori);
EdgePts_P(indreplace,:)=EdgePts(indreplace,:);

if iteration>1
    for i=1:iteration-1
        dis_cz_ori = sqrt(sum((midp-EdgePts_P) .^ 2,2));
        
        EdgePts_P1= EdgePts_P(1,:);
        EdgePts_P2 = (EdgePts_P(1:end-2,:)+EdgePts_P(2:end-1,:)+EdgePts_P(3:end,:))/3;
        EdgePts_P3 = EdgePts_P(end,:);
        EdgePts_P_ir = [EdgePts_P1;EdgePts_P2;EdgePts_P3];
        
        dis_midp_ir = sqrt(sum((midp-EdgePts_P_ir) .^ 2,2));
        indreplace = find(dis_midp_ir>=dis_cz_ori);
        EdgePts_P(indreplace,:)=EdgePts_P_ir(indreplace,:);
        
    end
end

return