function  UniquePtEdges_new = find_UniquePtEdges_3and4(NodeCount,UniquePtEdges)
% This is function is used to fix those points connnected with muitiple points

nodes4 = find(NodeCount==4);
x =zeros(size(UniquePtEdges));
ind1 = find(UniquePtEdges==nodes4);
nodes3 = find(NodeCount==3);
ind = [];
for j=1:size(nodes3,2)
    ind=[ind;find(UniquePtEdges==nodes3(j))]
end
x([ind;ind1])=1;
inddele = find(sum(x,2)==2)
UniquePtEdges(inddele,:)=[];
UniquePtEdges_new = UniquePtEdges;

return