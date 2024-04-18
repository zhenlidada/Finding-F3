function [Xpoint,inx] = find_Xpoint(EdgePts,X)
%====================================================================================================================
% This function is used to find the end point that is along with head mesh arch whithin certain distance.
% [input]:
%       EdgePts: a series of points created by the intersection plane of three points and a head mesh
%                the first row of EdgePts is the coordinate of the starting
%                point.
%       X:       the distance from the starting point to the end point.
% [output]
%       Xpoint:  the coordinate of ending point.
%       inx:     the index of the point beforehand the endign point.
% ====================================================================================================================
% Zhen Li, 2023, Maastricht university
% zhen.li@maastrichtuniversity.nl
% zhen.li.dada@gmail.com
% =====================================================================================================================
%%
cum_sum  = cumsum(sqrt(sum(diff(EdgePts).^2,2)));
inx = find(cum_sum <= X);
if ~isempty(inx) & ~isempty(max(inx)<size(EdgePts,1)-1|inx<=2)
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
else
    inx=0;
    dis_last = cum_sum(inx+1);
    last_more = X ;
    y = (EdgePts(inx+2,:)-EdgePts(inx+1,:))*(last_more/dis_last);
    x = EdgePts(inx+1,:);
    Xpoint = [x(1)+y(1),x(2)+y(2),x(3)+y(3)];
end
end
