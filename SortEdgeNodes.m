%% Paolo Giacometti
% Function takes the list of node coordinates in EdgeNodes and the matrix
% of node indices that compose it edge in Edges to find the unique set of
% points and sort them starting with the point closest to fiducial 1 and
% arranged in the direction towards fiducial 2. First, it finds point
% closest to fiducial 1 and direction of sorting by finding the point next
% to the one closest to fiducial 1 that is closest to fiducial 2. From
% there, the function finds the list of indices of repeating nodes. It
% iterates through all the edges between two points and organizes the
% points by finding the pairings of each point. Once the edges are sorted
% in positional order, the first column is selected as the indices of the
% nodes in order. From that set the list of sorted nodes is computed. If
% there are multiple loops of points in the set of edges, the function
% reruns the sorting on only the loop that contains the most amount of
% points. (This assumes that the nodes of interest are those in the largest
% loop)
% Inputs: Edges [e1 e2] list of indices mapping to two node points.
%         EdgeNodes [x,y,z] coordinates of all points.
%         Fiducials [x,y,z] coordinates of reference points
% Output: SortNodes [x,y,z] coordinates of unique and sorted points.
function [SortNodes] = SortEdgeNodes(Edges,EdgeNodes,Fiducials)

% Find Unique set of nodes
[UniqueNodes,IA,IC] = unique(EdgeNodes,'rows','stable');

% Regenerate set of edges that point to the unique set of nodes
UniquePtEdges = unique(IC(Edges),'rows','stable');

% Count the instances each node is in edges
NodeCount = hist(UniquePtEdges(:),1:size(UniqueNodes,1));

% Remove extra nodes and edges that aren't part of the loop
%     if ~isempty(find(NodeCount<2))

if ~isempty(find(NodeCount==0))
    warning('SortEdgeNodes warning:: Unconnected edge found. The mesh may have unconnected elements.')
    warning('SortEdgeNodes warning: NodeCount==0.')
    % Get index of extra nodes
    ExtraNodeidx = find(NodeCount==0);
    
    for i=1:numel(ExtraNodeidx)
        % Remove edges with extra nodes from list
        ExtraEdgeRem = UniquePtEdges(~sum(UniquePtEdges==ExtraNodeidx(i),2),:);
        
        % Change edges to reflect removed node
        ExtraEdgeRem(ExtraEdgeRem>ExtraNodeidx(i)) = ExtraEdgeRem(ExtraEdgeRem>ExtraNodeidx(i))-1;
        
        % Set new edges list as current
        UniquePtEdges = ExtraEdgeRem;
        
        % Remove node from list
        UniqueNodes = UniqueNodes([1:ExtraNodeidx(i)-1 ExtraNodeidx(i)+1:size(UniqueNodes,1)],:);
    end
end
%% fix for UniqueNodes
if ~isempty(find(NodeCount==1))
    warning('SortEdgeNodes warning:: NodeCount==1.')
    % Get index of extra nodes
    ExtraNodeidx = find(NodeCount==1);
    ind_rm = find(UniquePtEdges==ExtraNodeidx);
    for i=1:size(ind_rm)
        if ind_rm(i)>size(UniquePtEdges,1)
            removeind(i) = ind_rm(i)-size(UniquePtEdges,1);
        else
            removeind(i) = ind_rm(i);
        end
    end
    UniquePtEdges(removeind,:)=[];
end

nodes = find(NodeCount==4);
num = size(nodes,2);
node3 = find(NodeCount==3);

if ~isempty(node3)
    size3 = size(find(UniquePtEdges==node3),1);
end

if num>=1 & isempty(find(NodeCount==3))
    for i=1:num
        F1_4 = UniquePtEdges(find(sum(UniquePtEdges==nodes(i),2)),:);
        F1_4 = F1_4(F1_4~=nodes(i));
        inum(i,1) = size(F1_4,1);
    end
    ExtraNodeidx = find(NodeCount==4);
    ExtraNodeidx = ExtraNodeidx(find(inum==4));
    UniquePtEdges = find_UniquePtEdges_4(ExtraNodeidx,UniquePtEdges,UniqueNodes);
    
elseif ~isempty(node3)
    
    if size3==2
        for i=1:num
            F1_4 = UniquePtEdges(find(sum(UniquePtEdges==nodes(i),2)),:);
            F1_4 = F1_4(F1_4~=nodes(i));
            inum(i,1) = size(F1_4,1);
        end
        if num~=0
            ExtraNodeidx = find(NodeCount==4);
            ExtraNodeidx = ExtraNodeidx(find(inum==4));
            UniquePtEdges = find_UniquePtEdges_4(ExtraNodeidx,UniquePtEdges,UniqueNodes);
        end
    elseif size3>2&num==0
        UniquePtEdges = find_UniquePtEdges_3(NodeCount,UniquePtEdges,UniqueNodes);
    else
        UniquePtEdges = find_UniquePtEdges_3and4(NodeCount,UniquePtEdges);
    end
elseif ~isempty(find(NodeCount==3))&~isempty(find(NodeCount==4))
    UniquePtEdges = find_UniquePtEdges_3and4(NodeCount,UniquePtEdges);
end
%%
% Find point in nodes closest to first fiducial
DistP1 = sqrt(sum((UniqueNodes-repmat(Fiducials(1,:),size(UniqueNodes,1),1)).^2,2));
[DistP11,ind] = sort(DistP1);
i=1;
while isempty(find(UniquePtEdges==ind(i),1))
    i=i+1;
end
mpt=DistP1(ind(i));

% mpt = min(DistP1);

% Find index in nodes of that point
mpidx = find(DistP1==mpt);
% mpidx = min(find(DistP1==mpt));

% Find edge pairs to that point
MinEd = UniquePtEdges(sum(UniquePtEdges==mpidx,2)>0,:);

% Find edge index to those pairs
MinEdidx = find(sum(UniquePtEdges==mpidx,2)>0);

% Find edge point of those pairs
MinPairPt = UniqueNodes(MinEd(MinEd~=mpidx),:);
% x = MinEd(MinEd~=mpidx);
% MinPairPt = UniqueNodes(min(x),:);

% Find points in nodes closest to second fiducial
DistP2 = sqrt(sum((MinPairPt-[Fiducials(2,:);Fiducials(2,:)]).^2,2));
[mpt midx] = min(DistP2);
% [mpt midx] = max(DistP2);

% Initialize new edges matrix
SortEdges = zeros(size(UniquePtEdges));

% Replace index of duplicate nodes in Edges matrix with the index of
% the first occurrence.
% Select current edges connected to starting point

CurEdidx = MinEdidx(midx);
CurEd = UniquePtEdges(CurEdidx,:);
CurNode=CurEd(CurEd~=mpidx); % Select node connected to starting point

% Make new sorted Edges list starting with first edge
SortEdges(1,:) = [mpidx CurNode];

for i=2:size(UniquePtEdges,1)
    % Find edge indices connected to that pair
    NextEdgidx = find(sum(UniquePtEdges==CurNode,2)>0);
    
    % Select next edge (the one that isnt the current edge)
    NextEdge = UniquePtEdges(NextEdgidx(NextEdgidx~=CurEdidx),:);
    
    % Find edge index of point paired to the current point
    NewEdgPtidx = NextEdge(NextEdge~=CurNode);
    % Make new sorted Edges list
    SortEdges(i,:) = [CurNode NewEdgPtidx];
    % Select new edge point to find duplicates for
    CurNode = NewEdgPtidx;
    
    % Make New edge index current edge index for next iteration
    CurEdidx = NextEdgidx(NextEdgidx~=CurEdidx);
    
end

% Check to make sure only one loop was present in the edge
if size(unique(SortEdges(:,1)),1)~=size(SortEdges,1)
    
    % Return warning regarding extra loops.
    warning('SortEdgeNodes warning:: Multiple loops found. The mesh was sliced at a location that produced multiple edges.')
    % Find the amount of nodes on each loop
    if size(find(SortEdges(:,1)==SortEdges(1)),1)>2
        SortEdges = SortEdges(1:size(unique(SortEdges(:,1)),1),:);
        SortNodes = UniqueNodes(SortEdges(:,1),:);
    else
        
        LoopSz1 = size(unique(SortEdges(:,1)),1);
        LoopSz2 = size(SortEdges,1)-size(unique(SortEdges(:,1)),1);
        
        if LoopSz1>=LoopSz2  %LoopSz1>LoopSz2
            
            % Return only the edges belonging to that loop
            SortEdges = SortEdges(1:LoopSz1,:);
            
            % Pick nodes in sorted order from edges
            SortNodes = UniqueNodes(SortEdges(:,1),:);
            
        else
            % Get logical of points in larger loop
            LoopNodesLog = true(size(SortEdges,1),1);
            LoopNodesLog(unique(SortEdges(:,1))) = false;
            
            % Select Loop Nodes
            LoopNodes = UniqueNodes(LoopNodesLog,:);
            
            % Make list of new indices
            NewLoopNodeidx = (1:size(LoopNodes))';
            
            % Make set of indices with the size of original set
            NewLoopNodeOrigidx = zeros(size(UniqueNodes,1),1);
            
            % Populate set with new indices
            NewLoopNodeOrigidx(LoopNodesLog) = NewLoopNodeidx;
            
            % Get logical edges from larger loop
            LoopEdgLog = LoopNodesLog(UniquePtEdges(:,1));
            
            % Get original Edges in loop
            LoopEdgesOrig = UniquePtEdges(LoopEdgLog,:);
            
            % Make new set of edges for nodes in loop
            LoopEdges = NewLoopNodeOrigidx(LoopEdgesOrig);
            
            % Sort nodes within the loop
            SortNodes = SortEdgeNodes(LoopEdges,LoopNodes,Fiducials);
            
        end
    end
else
    
    % Pick nodes in sorted order from edges
    SortNodes = UniqueNodes(SortEdges(:,1),:);
    
end

return