function [a ,b, out] = maxx(A)
%% Returns the index and value of the maximum element of an array
    out = -inf;
    for i=1:size(A,1)
        for j=1:size(A,2)
            if(A(i,j)>out)
                out = A(i,j);
                a = i;
                b = j; 
            end
        end
    end
end