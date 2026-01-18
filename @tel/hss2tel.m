function T = hss2tel(H)
%Convert an @hss matrix in @tel format using the standard telescopic
%decomposition*
%INPUT: H matrix in @hss format
%OUTPUT: T matrix in @tel format
%*"Computing functions of symmetric hierarchically semiseparable matrices"
if ~isa(H,'hss') 
    error('Input must be an hm_toolbox @hss object.');
end
d = depth(H);
if d == 0 % H has no @hss structure
    T = H.D;
else
    T = tel;  
    T.top = 0;
    T = enter_level(T, H, 0, d);
end
end
%------------------------------------------------------------------------
function T = enter_level(T,H,c_depth,m_depth)
    %Builds level c_depth of T @tel with data in H @hss at level
    %m_depth-c_depth
    %INPUT: T @tel matrix, H @hss matrix, c_depth, m_depth the depth of H
    %OUTPUT: T @tel matrix, the conversion of H @hss 
    T = fill_level(T,H,1,c_depth,m_depth);
    if  c_depth + 1 == m_depth
        T.top = 1;
        T.M = [zeros(size(H.B12,1), size(H.B21,2)), H.B12; H.B21, zeros(size(H.B21,1),size(H.B12,2))];
    else  
        T.top = 0;
        T.M = tel;
        T.M = enter_level(T.M,H,c_depth+1,m_depth);
    end
end
%------------------------------------------------------------------------
function [T,index] = fill_level(T,H,index,c_depth,t_depth)
%Fills level c_depth U,V,D cells of T @tel with data from H @hss
%INPUT: T @tel, H @hss, index numbers items at level t_depth in H, c_depth depth
%of current level in @hss H navigation, t_depth depth of target level.
%OUTPUT: T @tel with level c_depth filled, index 
if c_depth < t_depth
    [T,index] = fill_level(T,H.A11,index, c_depth+1, t_depth);
    [T,index] = fill_level(T,H.A22,index, c_depth+1, t_depth);
else
    if H.leafnode == 1
        T.D{index} = H.D;
        T.U{index} = H.U;
        T.V{index} = H.V;
    else
    T.D{index} = [zeros(size(H.B12,1), size(H.B12,1)), H.B12; H.B21, zeros(size(H.B12,2),size(H.B12,2))];
    T.U{index} = [H.Rl; H.Rr];
    T.V{index} = [H.Wl; H.Wr];  
    end
    index = index+1;
end
end
% %--------------------------------------------------------------------------
% % function tot = total_size(C, dim)
% % % Computes cumulative size of a cell array of matrices.
% % %INPUT: cell array of matrices C.
% % %OUTPUT: tot=[tot_r,tot_c] cumulative rows and columns sizes.
% % tot = 0;
% % for i = 1:length(C)
% %     current = C{i};
% %     tot = tot + size(current,dim);
% % end
% % end
% %--------------------------------------------------------------------------
% function T = hss2tel(H)
%     d = depth(H);
%     if d == 0 
%         T = H.D; 
%         return; 
%     end
%     % 1. Raccogli tutti i nodi per livello in un unico passaggio
%     levels = cell(1, d + 1);
%     levels = collect_nodes(H, 0, levels); 
%     % 2. Costruisci la struttura tel iterando sui livelli raccolti
%     T = build_tel_recursive(levels, 0, d);
% end
% 
% function T = build_tel_recursive(levels, curr_d, max_d)
%     T = tel;
%     nodes = levels{max_d - curr_d + 1}; % Naviga i livelli HSS
%     num_nodes = length(nodes);
%     T.D = cell(1, num_nodes); 
%     T.U = cell(1, num_nodes); 
%     T.V = cell(1, num_nodes);
%     for i = 1:num_nodes
%         h = nodes{i};
%         if h.leafnode
%             T.D{i} = h.D; T.U{i} = h.U; T.V{i} = h.V;
%         else
%             %Proposizione 3.2: blocchi core fuori diagonale
%             T.D{i} = [zeros(size(h.B12,1)), h.B12; h.B21, zeros(size(h.B21,1))];
%             T.U{i} = [h.Rl; h.Rr]; T.V{i} = [h.Wl; h.Wr];
%         end
%     end
%     if curr_d < max_d
%         T.top = 0;
%         T.M = build_tel_recursive(levels, curr_d + 1, max_d);
%     else
%         T.top = 1;
%         % Il blocco più interno (M denso)
%         T.M = zeros(0); % O il core block della radice
%     end
% end

function d = depth(H)
    % Returns the depth of the perfect binary tree associated to H @hss
    if H.leafnode
        d = 0;     %root has depth 0        
    else
        d = 1 + depth(H.A11);  
    end
end