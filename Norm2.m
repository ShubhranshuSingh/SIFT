function output = Norm2(M)
%% Returns 2-Norm of a vector
    rows = size(M,1);
    cols = size(M,2);
    output = 0;
   for i=1:rows
        for j=1:cols
            output = output + M(i,j)^2;
        end
   end
   output = sqrt(output);
end    