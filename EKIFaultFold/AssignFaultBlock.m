function block = AssignFaultBlock(x,y,z,faults,fault_blocks_relationships)
%Assign fault block numbers to points.
%x, y, and z should be column vectors.
%faults should be an array of fault_class objects.
%fault_blocks_relationships should be a 2D array as in the options file.

block = zeros(size(x));
for j = 1:size(fault_blocks_relationships,1) %Look through fault blocks.
    active = true(size(x));
    for k = 1:size(fault_blocks_relationships,2) %Loop through faults.
        if fault_blocks_relationships(j,k) ~= 0
            bed_uv = faults(k).M\[x,y,z,ones(size(x))]';
            f = InterpolateSurface(faults(k),bed_uv(1,:),bed_uv(2,:));
            hw = bed_uv(3,:)>=f;
            fw = bed_uv(3,:)<f;
            if fault_blocks_relationships(j,k) == -1 %Fault block is in the footwall.
                active(hw) = false;
            elseif fault_blocks_relationships(j,k) == 1 %Fault block is in the hanging wall.
                active(fw) = false;
            end
        end
    end
%     if any(block(active)~=0)
%         warning([num2str(sum(block(active)~=0)),' points may be assigned to multiple fault blocks. Using first one.'])
%     end
    block(active & block==0) = j;
end
if any(block==0)
    warning(['Fault blocks were not assigned to ',num2str(sum(block==0)),' points.'])
end

end