 function [Q] = randPDsymmMatix(n)
 % This function generates random positive definite symmetric matrices of size n. 
 % By Pablo Nanez (Ñañez)
 Q = randn(n);
 Q = Q*Q';
 Q = symmetry(Q,1);
 fprintf('The eigenvalues\n')
 eig(Q)
 end

 function [A] = symmetry(A,tipo)
 % Funcion para hacer la matriz simetrica de A
 % Tipo: 1: matriz simetrica con signos
 % Tipo: -1: matriz simetrica con signos contrarios
 % e.g., D = simetria(D,1);
 [f c] = size(A);
 for iff = 1:f
     for ic = 1:c
         if isempty(A(iff,ic)) || A(iff,ic)==0
         else
             A(ic,iff) = tipo*A(iff,ic);
         end
     end
 end
 end 