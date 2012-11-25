function [MemberIndex,GroupCenters,GroupSizes,NumGroups,result]=...
    groupScalar(x,radius)
%GROUPSCALAR group scalars according to their similarity.
% radius: group radius (half width in values of scalars of the same group)
%
% returns:
% MemberIndex: MaxNumMembers x NumGroups, indices of group members
% GroupCenters: group centers (arithmetic average)
% GroupSizes:  1 x NumGroups, Number of members of each group
% NumGroups:   number of groups

d=radius*2;

MaxNumGroups=ceil((max(x(:))-min(x(:)))/d);

[x,ix]=sort(x(:));
L=length(x(:));

MemberIndex=cell(1,MaxNumGroups);
GroupCenters=zeros(1,MaxNumGroups);
GroupSizes=zeros(1,MaxNumGroups);
GroupBots=zeros(1,MaxNumGroups);
GroupSums=zeros(1,MaxNumGroups);

g=1;
GroupBots(1)=x(1);
GroupSums(1)=x(1);
MemberIndex{1}=zeros(L,1);
MemberIndex{1}(1)=ix(1);
GroupSizes(1)=1;
for k=2:L
    if x(k)-GroupBots(g)<2*radius
        GroupSums(g)=GroupSums(g)+x(k);
        GroupSizes(g)=GroupSizes(g)+1;
        MemberIndex{g}(GroupSizes(g))=ix(k);
    else
        g=g+1;
        GroupBots(g)=x(k);
        GroupSums(g)=x(k);
        GroupSizes(g)=1;
        MemberIndex{g}=zeros(L-k+1,1);
        MemberIndex{g}(GroupSizes(g))=ix(k);
    end
end
NumGroups=g;

GroupCenters(1:NumGroups)=GroupSums(1:NumGroups)./GroupSizes(1:NumGroups);

result={MemberIndex,GroupCenters,GroupSizes,NumGroups,x,radius};
return
