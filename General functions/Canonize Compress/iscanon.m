function canon = iscanon(mps)
%%
canon = 0;
tolerance = 1E-12;
isleft = 0;
isright = 0;

N = length(mps);
% is left?

for i = 1:N
    s = size(mps{i});
    M = mps{i};
    test = contract(M,[1,3],conj(M),[1,3]);
    if ~approx(max(max(test-eye(s(2)))),0,tolerance)
        break
    end
    if i == N
        isleft = 1;
    end
end
% is right?

if ~isleft
    for i = 1:N
        s = size(mps{i});
        M = mps{i};
        test = contract(M,[2,3],conj(M),[2,3]);
        if ~approx(max(max(test-eye(s(1)))),0,tolerance)
            break
        end
        if i == N
            isright = 1;
        end
    end
end

if isleft
    canon = 1;
end
if isright
    canon = -1;
end
end

