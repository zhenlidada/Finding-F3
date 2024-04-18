function UniquePtEdges_new = find_UniquePtEdges_one_4(NodeCount,UniquePtEdges)
% This is function is used to fix those points connnected with muitiple points

ExtraNodeidx = find(NodeCount==4);
F1_4 = UniquePtEdges(find(sum(UniquePtEdges==ExtraNodeidx,2)),:);
F1_4 = F1_4(F1_4~=ExtraNodeidx);

for i =1:4
    nextone = F1_4(i);
    ran1_all = ExtraNodeidx;
    ran1_all = [ran1_all,nextone(1)];
    k=2;
    
    while ~isempty(nextone)
        ran1_all(k,1) = nextone(1);
        nextone = UniquePtEdges(find(sum(UniquePtEdges==nextone,2)),:);
        nextone = nextone(nextone~=ran1_all(k,1)& ~ismember(nextone,ran1_all));
        if ~isempty(nextone)
            ran1_all(k,2)=nextone;
            ran1_all(k+1,1)=ExtraNodeidx;
            ran1_all(k+1,2)=ExtraNodeidx;
            k=k+1;
        else
            ran1_all(k,2)=ExtraNodeidx;
            break
            
        end
    end
    
    Loop_all{i,1} = ran1_all;
    Loop_all{i,2} = size(ran1_all,1);
    
end

[s1,s2,s3,s4] = (Loop_all{:,2});
ind_4 = [s1,s2,s3,s4];
if  max(ind_4)< size(UniquePtEdges,1)/2
    indlast_end=[];
    for j=1:size(Loop_all{1,1},1)
        indlast_end = [indlast_end;find(sum(UniquePtEdges==Loop_all{1,1}(j,:),2)==2|sum(UniquePtEdges==flip(Loop_all{1,1}(j,:),2),2)==2,1)];
    end
    UniquePtEdges(indlast_end,:)=[];
    UniquePtEdges_new = UniquePtEdges;
    
    if size(unique(UniquePtEdges_new),1)~=size(UniquePtEdges_new,1)
        UniquePtEdges_new = add_edge(UniquePtEdges_new);
    end
else
    ind_loop = find(ind_4==max(ind_4),1);
    UniquePtEdges_new = Loop_all{ind_loop};
end

return

