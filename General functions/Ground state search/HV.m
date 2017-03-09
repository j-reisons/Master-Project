function [handle] = HV(Right,Left,Operator)
%HV Returns function handle for use with eigenvalue solver

    s_O = size(Operator);
    s_R = size(Right);
    s_L = size(Left);
    
    d = s_O(3);
    D_L = s_L(1);
    D_R = s_R(1);

    function [Hx] = Hvfun(x)
        M = reshape(x,[D_L,D_R,d]);
        
        if Right == 1
            R_O = reshape(Operator,[1,1,s_O(1),s_O(3),s_O(4)]);
        else
        R_O = contract(Right,2,Operator,2);
        end
        
        R_O_M = contract(R_O,[1,5],M,[2,3]);
        
        if Left == 1
            s = size(R_O_M);
            R_O_M_L = reshape(R_O_M,[1,s(1),s(3)]);
            
        else
            R_O_M_L = contract(R_O_M,[4,2],Left,[1,2]);
            R_O_M_L = permute(R_O_M_L,[3,1,2]);
        end
        
        Hx = reshape(R_O_M_L,[D_L*D_R*d,1]);
    end

handle = @Hvfun;
end

