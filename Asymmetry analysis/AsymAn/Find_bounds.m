function [ p_min, p_max ] = Find_bounds(p, norm_K2 )
%FIND_BOUNDS Summary of this function goes here
%   Detailed explanation goes here

% Find i_opt
[~,i_opt] = min(norm_K2);

% Find the first value that goes over 1 in +ve direction
norm_Ki = 0;
i = 0;
Out = false;

while norm_Ki <= 1 && Out == false
    i = i+1;
    if i_opt+i <= size(norm_K2,1)
        norm_Ki = norm_K2(i_opt+i);
    else
        Out = true;
    end
end

if Out == false
    p_max = p(i_opt+i);
else
    p_max = NaN;
end


% Find the first value that goes over 1 in -ve direction
norm_Ki = 0;
i = 0;
Out = false;

while norm_Ki <= 1 && Out == false
    i = i-1;
    if i_opt+i > 0
        norm_Ki = norm_K2(i_opt+i);
    else
        Out = true;
    end
end



if Out == false
    p_min = p(i_opt+i);
else
    p_min = NaN;
end

end

