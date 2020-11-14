function [ SmallObjects, BigObjects ] = SieveObject(MinSize, MaxSize, ParamObject )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n_object = size(ParamObject,1);
SmallObjects = zeros(1,n_object);
BigObjects = zeros(1,n_object);


for i = 1:n_object
    if ParamObject(i).Area > MaxSize || ParamObject(i).Area < MinSize
        SmallObjects(i) = i;
    else
        BigObjects(i) = i;
    end
end

SmallObjects(SmallObjects == 0) = [];
BigObjects(BigObjects == 0) = [];


idx = find([ParamObject.Area] > MinSize & [ParamObject.Area] < MaxSize & [ParamObject.Eccentricity] < 1)
% z = ismember(L,idx);

end

