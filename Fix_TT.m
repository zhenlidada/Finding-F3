function [arcLength,PtSort] = Fix_TT(arcLength,PtSort)
% fix Arch path from left trgus and right tragus during getting Cz
% by remove the points inside of the acoustic meatus

inclu1 = find(diff(PtSort(1:round(size(PtSort,1)/4),3))<0);
if ~isempty(inclu1)
    if ~isempty(find(inclu1~=1))
        if inclu1(1)==1
            inclu1 = inclu1(2);
        else
            inclu1 = inclu1(1);
        end
        exclu1 = inclu1+1;
        inclu1_co = PtSort(inclu1,:);
        exerror = 5;
        diffex = sqrt(sum((inclu1_co-PtSort(exclu1+exerror:end,:)).^2,2));
        inclu2 = find(diffex==min(diffex))+inclu1;
        PtSort2 = [PtSort(1:inclu1,:);PtSort(inclu2+exerror:end,:)];
        if size(PtSort2,1) ~= size(PtSort,1)
            arcLength =cumsum(sqrt(sum((diff(PtSort2).^2),2)));
            arcLength = arcLength(end);
            PtSort = PtSort2;
        end
    end
end

PtSort = flip(PtSort,1);
inclu1 = find(diff(PtSort(1:round(size(PtSort,1)/4),3))<0);
if ~isempty(inclu1)
    if ~isempty(find(inclu1~=1))
        if inclu1(1)==1
            inclu1 = inclu1(2);
        else
            inclu1 = inclu1(1);
        end
        exclu1 = inclu1+1;
        inclu1_co = PtSort(inclu1,:);
        exerror = 5;
        diffex = sqrt(sum((inclu1_co-PtSort(exclu1+exerror:end,:)).^2,2));
        inclu2 = find(diffex==min(diffex))+inclu1;
        PtSort2 = [PtSort(1:inclu1,:);PtSort(inclu2+exerror:end,:)];
        if size(PtSort2,1) ~= size(PtSort,1)
            arcLength =cumsum(sqrt(sum((diff(PtSort2).^2),2)));
            arcLength = arcLength(end);
            PtSort = PtSort2;
        end
    end
end

PtSort = flip(PtSort,1);

return