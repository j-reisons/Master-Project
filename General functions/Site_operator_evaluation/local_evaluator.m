function [evaluations] = local_evaluator(mps,canonization,varargin)
%LOCAL_EVALUATOR Efficiently evaluates n_site operators on canonized states
Operators = varargin;
N = length(mps);
s = size(mps{1});
d = s(3);

N_ops = length(Operators);
evaluations = cell(1,N_ops);

op_sites = zeros(N_ops);
for i = 1:N_ops
    op_sites(i) = length(Operators{i});
    evaluations{i} = zeros(1,N - op_sites(i) + 1);
end

switch canonization
    case -1
        
        for j = 1:N_ops
        mpo = Operators{j};
        eval = contract(conj(mps{1}),3,mpo{1},4);
        eval = contract(eval,5,mps{1},3);
        eval = squeeze(eval);
        evaluations{j}(1) = trace(eval);
        end
        
        L = squeeze(contract(mps{1},3,conj(mps{1}),3));
        
        for i = 2:1:N
            L = contract(L,1,mps{i},1);
            L = contract(L,1,conj(mps{i}),1);
            for j = find(op_sites < N-i)
                eval = L;
                eval = contract(eval,[4,2],mpo{i},[3,4]);
                eval = squeeze(eval);

                if op_sites(j) > 1 % For n > 1 operators
                    for k = 1 : op_sites(j) - 1
                        eval = contract(eval,1,mps{i+k},1);
                        eval = contract(eval,1,conj(mps{i+k}),1);
                        eval = contract(eval,[4,2],mpo{i+k},[3,4]);
                    end
                end
                evaluations{j}(i) = trace(eval);
                
            end
            
        L = contract(L,[4,2],eye(d),[1,2]); %Lazy internal contraction
        end
        
        
    case 1
end

end

