 function [x,y,z] = FilterHWFW(x,y,z,active,FaultBlock,faults,fault_blocks_relationships)
%Filter out points that are in the wrong fault block.

%This script is not actually used by anything under RunEKI. However, I have
%used it in some of my plotting functions.

nfaults = length(faults);
mask = true(size(active)); %Initialize as all true and figure out which ones should be active.
for j = 1:nfaults
    %Transform the bed points to the rotated coordinate system.
    bed_uv = faults(j).M\[x,y,z,ones(size(x))]'; 
    f = InterpolateSurface(faults(j),bed_uv(1,:),bed_uv(2,:));
    hw = bed_uv(3,:)>=f;
    fw = bed_uv(3,:)<=f;
    for k = 1:size(fault_blocks_relationships,1)
        block = k-1; %Fault blocks are numbered starting from 0.
        if fault_blocks_relationships(k,j) == 1
            mask((FaultBlock==block & fw')) = false;
        elseif fault_blocks_relationships(k,j) == -1
            mask(FaultBlock==block & hw') = false;
        end
    end
end
[x,y,z] = deal(x(mask),y(mask),z(mask));
end