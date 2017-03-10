function [evaluations] = Canon_evaluator(mps,canon,varargin)
% Evaluates n_site operators efficiently by mixed canonization
% As implemented, the contraction is not very efficient for large

Operators = varargin;
N = length(mps);
s = size(mps{1});
d = s(3);

N_ops = length(Operators);
evaluations = cell(1,N_ops);

op_sites = zeros(1,N_ops);

for i = 1:N_ops
    op_sites(i) = length(Operators{i});
    evaluations{i} = zeros(1,N - op_sites(i) + 1);
end

switch canon
    case -1
        for i = 1:N
            for j = find(op_sites <= N - i + 1)
                sites = op_sites(j);
                mpo = Operators{j};
                
                mini_mps = cell(1,sites);
                
                for k = 1:sites
                    mini_mps{k} = mps{i+(k-1)};
                end
                
                mini_mpo_mps = apply(mpo,mini_mps);
                
                eval = contract(mini_mpo_mps{1},[1,3],conj(mini_mps{1}),[1,3]);
                for k = 2:sites
                    eval = contract(eval,1,mini_mpo_mps{k},1);
                    eval = contract(eval,[1,3],conj(mini_mps{k}),[1,3]);
                end
                eval = squeeze(eval);
                evaluations{j}(i) = trace(eval);      
            end
            mps = L_can(mps,i);
        end
        
    case 1
        
        for i = N:-1:1
            for j = find(op_sites <= N - i + 1)
                sites = op_sites(j);
                mpo = Operators{j};
                mini_mps = cell(1,sites);

                for k = 1:sites
                    mini_mps{k} = mps{i+(k-1)};
                end
                
                mini_mpo_mps = apply(mpo,mini_mps);
                eval = contract(mini_mpo_mps{1},[1,3],conj(mini_mps{1}),[1,3]);
                for k = 2:sites
                    eval = contract(eval,1,mini_mpo_mps{k},1);
                    eval = contract(eval,[1,3],conj(mini_mps{k}),[1,3]);
                end
                eval = squeeze(eval);
                evaluations{j}(i) = trace(eval);      
            end
            mps = R_can(mps,i);
        end
        
        
end

end

