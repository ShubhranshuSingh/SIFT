function out = padreflect(I)
%% Retuns a padded image using reflection padding
    out = zeros(size(I,1)+2,size(I,2)+2);
    
    out(2:end-1,2:end-1) = I;
    
    out(1,:) = [I(1,1) I(2,:) I(1,end)];
    out(2:end,end) = [I(:,end-1); I(end,end)];
    out(end,1:end-1) = [I(end,1) I(end-1,:)];
    out(2:end-1,1) = I(:,2);
end