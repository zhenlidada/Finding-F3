function UniquePtEdges_new = add_edge(UniqueNodes,UniquePtEdges_new)
% This function is used to add one edge that missing between to point.

all_new = sort(UniquePtEdges_new(:));
q1 = 1:2:size(all_new,1)-1;
q2 = 2:2:size(all_new,1);
diff = all_new(q1)-all_new(q2);
ind = find(diff~=0);

if isempty(ind)
    NodeCount = hist(UniquePtEdges_new(:),1:size(UniqueNodes,1));
    if find(NodeCount==4)
        sizep = size(find(NodeCount==4),1);
        if sizep == 1
            UniquePtEdges_new = find_UniquePtEdges_4(NodeCount,UniquePtEdges_new);
        end
    end
else
    
    ind = ind(end)*2-1;
    all_new_egde = [];
    all_new_egde = [all_new_egde;all_new(ind)];
    
    while ~isempty(ind)
        all_new(ind)=[];
        q1 = 1:2:size(all_new,1)-1;
        q2 = 2:2:size(all_new,1);
        diff = all_new(q1)-all_new(q2);
        ind = find(diff~=0);
        if isempty(ind)
            break
        end
        ind = ind(1)*2-1;
        all_new_egde = [all_new_egde;all_new(ind)];
    end
    
    row1=1:2:size(all_new_egde,1)-1;
    row2=2:2:size(all_new_egde,1);
    all_new_egde_last(1:size(all_new_egde,1)/2,1) = all_new_egde(row1);
    all_new_egde_last(1:size(all_new_egde,1)/2,2) = all_new_egde(row2);
    if size(find(UniquePtEdges_new==all_new_egde(row1)))==1
        UniquePtEdges_new=[UniquePtEdges_new;all_new_egde_last];
    else
        ind_rm = find(sum(UniquePtEdges_new==all_new_egde_last,2)==2);
        UniquePtEdges_new(ind_rm,:)=[];
    end
    
    NodeCount = hist(UniquePtEdges_new(:),1:size(UniqueNodes,1));
    if find(NodeCount==4)
        sizep = size(find(NodeCount==4),1);
        if sizep == 1
            UniquePtEdges_new = find_UniquePtEdges_one_4(NodeCount,UniquePtEdges_new);
        end
    end
end

return