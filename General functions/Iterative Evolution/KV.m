function [handle] = KV(Right,Left)
%HV Returns function handle for zero-site effective Hamiltonian
    function [Kres] = Kvfun(res)    
        L_res = contract(res,1,Left,1);
        Kres = contract(L_res,[1,2],Right,[1,2]);
    end

handle = @Kvfun;
end

