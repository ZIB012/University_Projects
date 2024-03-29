%--------------------------------------------------------------------
% PURPOSE:
%
% This routine generates a structure containing the scalar-valued Lagrange
% shape functions (and their derivatives) for P^k discontinuous finite
% elements with k=1,2
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [basis]= basis_lagrange(fem)

switch fem
    case 'P1'
        nln=3;
        [c]= matrix_coeff_P1(nln);
        basis=struct('num',nln,...
            'n_edge',3,...
            'coeff',c,...
            'fbasis',{'c(1,1).*csi + c(2,1).*eta + c(3,1)',...
                      'c(1,2).*csi + c(2,2).*eta + c(3,2)',...
                      'c(1,3).*csi + c(2,3).*eta + c(3,3)'},...
            'Gbasis_1',{'c(1,1).*1 + 0     .*eta + 0',...
            'c(1,2).*1 + 0     .*eta + 0',...
            'c(1,3).*1 + 0     .*eta + 0'},...
            'Gbasis_2',{'0     .*csi + c(2,1).*1   + 0',...
            '0     .*csi + c(2,2).*1   + 0',...
            '0     .*csi + c(2,3).*1   + 0'});

    case 'P2'
        nln=6;
        [c]= matrix_coeff_P2(nln);
        basis=struct('num',nln,...
            'n_edge',6,...
            'coeff',c,...
            'fbasis',{'c(1,1).*csi.^2 + c(2,1).*eta.^2 + c(3,1).*csi.*eta + c(4,1).*csi + c(5,1).*eta  + c(6,1)',...
            'c(1,2).*csi.^2 + c(2,2).*eta.^2 + c(3,2).*csi.*eta + c(4,2).*csi + c(5,2).*eta  + c(6,2)',...
            'c(1,3).*csi.^2 + c(2,3).*eta.^2 + c(3,3).*csi.*eta + c(4,3).*csi + c(5,3).*eta  + c(6,3)',...
            'c(1,4).*csi.^2 + c(2,4).*eta.^2 + c(3,4).*csi.*eta + c(4,4).*csi + c(5,4).*eta  + c(6,4)',...
            'c(1,5).*csi.^2 + c(2,5).*eta.^2 + c(3,5).*csi.*eta + c(4,5).*csi + c(5,5).*eta  + c(6,5)',...
            'c(1,6).*csi.^2 + c(2,6).*eta.^2 + c(3,6).*csi.*eta + c(4,6).*csi + c(5,6).*eta  + c(6,6)'},...
            'Gbasis_1',{'c(1,1).*2.*csi + c(3,1).*eta + c(4,1)',...
            'c(1,2).*2.*csi + c(3,2).*eta + c(4,2)',...
            'c(1,3).*2.*csi + c(3,3).*eta + c(4,3)',...
            'c(1,4).*2.*csi + c(3,4).*eta + c(4,4)',...
            'c(1,5).*2.*csi + c(3,5).*eta + c(4,5)',...
            'c(1,6).*2.*csi + c(3,6).*eta + c(4,6)'},...
            'Gbasis_2',{'2.*c(2,1).*eta + c(3,1).*csi + c(5,1)',...
            '2.*c(2,2).*eta + c(3,2).*csi + c(5,2)',...
            '2.*c(2,3).*eta + c(3,3).*csi + c(5,3)',...
            '2.*c(2,4).*eta + c(3,4).*csi + c(5,4)',...
            '2.*c(2,5).*eta + c(3,5).*csi + c(5,5)',...
            '2.*c(2,6).*eta + c(3,6).*csi + c(5,6)'});

     case 'P3'
        nln=10;
        [c]= matrix_coeff_P3(nln);
        basis=struct('num',nln,...
            'n_edge',10,...
            'coeff',c,...
            'fbasis',{
            'c(1,1).*csi.^3 + c(2,1).*eta.^3 + c(3,1).*csi.^2.*eta + c(4,1).*csi.*eta.^2 + c(5,1).*csi.^2  + c(6,1).*eta.^2 + c(7,1).*eta.*csi + c(8,1).*csi + c(9,1).*eta +  c(10,1)',...
            'c(1,2).*csi.^3 + c(2,2).*eta.^3 + c(3,2).*csi.^2.*eta + c(4,2).*csi.*eta.^2 + c(5,2).*csi.^2  + c(6,2).*eta.^2 + c(7,2).*eta.*csi + c(8,2).*csi + c(9,2).*eta +  c(10,2)',...
            'c(1,3).*csi.^3 + c(2,3).*eta.^3 + c(3,3).*csi.^2.*eta + c(4,3).*csi.*eta.^2 + c(5,3).*csi.^2  + c(6,3).*eta.^2 + c(7,3).*eta.*csi + c(8,3).*csi + c(9,3).*eta +  c(10,3)',...
            'c(1,4).*csi.^3 + c(2,4).*eta.^3 + c(3,4).*csi.^2.*eta + c(4,4).*csi.*eta.^2 + c(5,4).*csi.^2  + c(6,4).*eta.^2 + c(7,4).*eta.*csi + c(8,4).*csi + c(9,4).*eta +  c(10,4)',...
            'c(1,5).*csi.^3 + c(2,5).*eta.^3 + c(3,5).*csi.^2.*eta + c(4,5).*csi.*eta.^2 + c(5,5).*csi.^2  + c(6,5).*eta.^2 + c(7,5).*eta.*csi + c(8,5).*csi + c(9,5).*eta +  c(10,5)',...
            'c(1,6).*csi.^3 + c(2,6).*eta.^3 + c(3,6).*csi.^2.*eta + c(4,6).*csi.*eta.^2 + c(5,6).*csi.^2  + c(6,6).*eta.^2 + c(7,6).*eta.*csi + c(8,6).*csi + c(9,6).*eta +  c(10,6)',...
            'c(1,7).*csi.^3 + c(2,7).*eta.^3 + c(3,7).*csi.^2.*eta + c(4,7).*csi.*eta.^2 + c(5,7).*csi.^2  + c(6,7).*eta.^2 + c(7,7).*eta.*csi + c(8,7).*csi + c(9,7).*eta +  c(10,7)',...
            'c(1,8).*csi.^3 + c(2,8).*eta.^3 + c(3,8).*csi.^2.*eta + c(4,8).*csi.*eta.^2 + c(5,8).*csi.^2  + c(6,8).*eta.^2 + c(7,8).*eta.*csi + c(8,8).*csi + c(9,8).*eta +  c(10,8)',...
            'c(1,9).*csi.^3 + c(2,9).*eta.^3 + c(3,9).*csi.^2.*eta + c(4,9).*csi.*eta.^2 + c(5,9).*csi.^2  + c(6,9).*eta.^2 + c(7,9).*eta.*csi + c(8,9).*csi + c(9,9).*eta +  c(10,9)',...
            'c(1,10).*csi.^3 + c(2,10).*eta.^3 + c(3,10).*csi.^2.*eta + c(4,10).*csi.*eta.^2 + c(5,10).*csi.^2  + c(6,10).*eta.^2 + c(7,10).*eta.*csi + c(8,10).*csi + c(9,10).*eta +  c(10,10)'},...
            'Gbasis_1',{
            'c(1,1).*3*csi.^2  + c(3,1).*2*csi.*eta + c(4,1).*eta.^2 + c(5,1).*2*csi   + c(7,1).*eta + c(8,1) ',...
            'c(1,2).*3*csi.^2  + c(3,2).*2*csi.*eta + c(4,2).*eta.^2 + c(5,2).*2*csi   + c(7,2).*eta + c(8,2) ',...
            'c(1,3).*3*csi.^2  + c(3,3).*2*csi.*eta + c(4,3).*eta.^2 + c(5,3).*2*csi   + c(7,3).*eta + c(8,3) ',...
            'c(1,4).*3*csi.^2  + c(3,4).*2*csi.*eta + c(4,4).*eta.^2 + c(5,4).*2*csi   + c(7,4).*eta + c(8,4) ',...
            'c(1,5).*3*csi.^2  + c(3,5).*2*csi.*eta + c(4,5).*eta.^2 + c(5,5).*2*csi   + c(7,5).*eta + c(8,5) ',...
            'c(1,6).*3*csi.^2  + c(3,6).*2*csi.*eta + c(4,6).*eta.^2 + c(5,6).*2*csi   + c(7,6).*eta + c(8,6) ',...
            'c(1,7).*3*csi.^2  + c(3,7).*2*csi.*eta + c(4,7).*eta.^2 + c(5,7).*2*csi   + c(7,7).*eta + c(8,7) ',...
            'c(1,8).*3*csi.^2  + c(3,8).*2*csi.*eta + c(4,8).*eta.^2 + c(5,8).*2*csi   + c(7,8).*eta + c(8,8) ',...
            'c(1,9).*3*csi.^2  + c(3,9).*2*csi.*eta + c(4,9).*eta.^2 + c(5,9).*2*csi   + c(7,9).*eta + c(8,9) ',...
            'c(1,10).*3*csi.^2 + c(3,10).*2*csi.*eta + c(4,10).*eta.^2 + c(5,10).*2*csi  +  c(7,10).*eta + c(8,10)'},...
            'Gbasis_2',{
            'c(2,1).*3*eta.^2  + c(3,1)*csi.^2 + c(4,1).*2*csi.*eta  + c(6,1).*2*eta   + c(7,1).*csi + c(9,1) ',...
            'c(2,2).*3*eta.^2  + c(3,2)*csi.^2 + c(4,2).*2*csi.*eta  + c(6,2).*2*eta   + c(7,2).*csi + c(9,2) ',...
            'c(2,3).*3*eta.^2  + c(3,3)*csi.^2 + c(4,3).*2*csi.*eta  + c(6,3).*2*eta   + c(7,3).*csi + c(9,3) ',...
            'c(2,4).*3*eta.^2  + c(3,4)*csi.^2 + c(4,4).*2*csi.*eta  + c(6,4).*2*eta   + c(7,4).*csi + c(9,4) ',...
            'c(2,5).*3*eta.^2  + c(3,5)*csi.^2 + c(4,5).*2*csi.*eta  + c(6,5).*2*eta   + c(7,5).*csi + c(9,5) ',...
            'c(2,6).*3*eta.^2  + c(3,6)*csi.^2 + c(4,6).*2*csi.*eta  + c(6,6).*2*eta   + c(7,6).*csi + c(9,6) ',...
            'c(2,7).*3*eta.^2  + c(3,7)*csi.^2 + c(4,7).*2*csi.*eta  + c(6,7).*2*eta   + c(7,7).*csi + c(9,7) ',...
            'c(2,8).*3*eta.^2  + c(3,8)*csi.^2 + c(4,8).*2*csi.*eta  + c(6,8).*2*eta   + c(7,8).*csi + c(9,8) ',...
            'c(2,9).*3*eta.^2  + c(3,9)*csi.^2 + c(4,9).*2*csi.*eta  + c(6,9).*2*eta   + c(7,9).*csi + c(9,9) ',...
            'c(2,10).*3*eta.^2 + c(3,10)*csi.^2 + c(4,10).*2*csi.*eta  + c(6,10).*2*eta  +  c(7,10).*csi + c(9,10)'});       
        
        
        
        
    otherwise
        error('fem can be only P1-P2-P3');
end


