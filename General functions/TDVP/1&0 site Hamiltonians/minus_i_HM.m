function [handle] = minus_i_HM(Right,Left,Operator)
%minus_i_HV Summary of this function goes here

    function [HM] = HMfun(M)
        
        if Right == 1
            s_O = size(Operator);
            R_O = reshape(Operator,[1,1,s_O(1),s_O(3),s_O(4)]);
        else
            R_O = contract(Right,2,Operator,2);
        end
        
        R_O_M = contract(R_O,[1,5],M,[2,3]);
        
        if Left == 1
            s = size(R_O_M);
            R_O_M_L = reshape(R_O_M,[1,s(1),s(3)]);
            
        else
            R_O_M_L = contract(Left,[1,2],R_O_M,[4,2]);
        end
        
        HM = -1i*R_O_M_L;
    end

handle = @HMfun;

end

