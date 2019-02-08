% Takes a list of atoms and their coordinates, and the basis set information 
% read in by basisread, and generates a full list of basis functions 
% for a given molecule.
%
% Input:
%   atoms         list of element numbers (array with K elements); 
%                 e.g. [6 8] for CO 
%   xyz_a0        K×3 array of Cartesian coordinates of nuclei, in bohr
%   basissetdef   basis set definitions, as read in by basisread
%   
% Output:
%   basis     M-element structure array, where each element contains:
%     basis(p).atom   element number (6 for C, 8 for O, etc.)     
%     basis(p).A      vector of Cartesian coordinates of the nucleus, in bohr
%     basis(p).a      vector of Cartesian exponents ([0 0 0] for s, [1 0 0] for px, etc.)
%     basis(p).alpha  array of radial exponents of primitives, in inverse bohr
%     basis(p).d      array of contraction coefficients
%     basis(p).N      array of normalization constants

function basis = buildbasis(atoms,xyz_a0,bdef)
    % Cartesian Exponents for each subshell.
    s    = [0 0 0];
    p_x  = [1 0 0];
    p_y  = [0 1 0];
    p_z  = [0 0 1];
    d_x2 = [2 0 0];
    d_xy = [1 1 0];
    d_xz = [1 0 1];
    d_y2 = [0 2 0];
    d_yz = [0 1 1];
    d_z2 = [0 0 2];
    
    %basis = struct;
    fac2 = @(n) prod(n:-2:1);

    for K=1:numel(atoms) % iterating over each atom in the atom list
        
        cart_exp = [];
        
        for i=1:numel(bdef{K}.shelltype) % iterating over shell types 
                                                % (i.e. 'S' then 'SP'..) 
            
            type = bdef{K}(i).shelltype;
            rad_exp = {};
            coeffs = {};
            q=1;
            ifSP = 1;
            
            if type == 'S'
                subshells = [s];
                X=0;
            elseif type == 'SP'
                subshells = [s; p_x; p_y; p_z];
                X=3;
                ifSP = 2;
            elseif type == 'P'
                subshells = [p_x; p_y; p_z];
                X=2;
            elseif type == 'D'
                subshells = [d_x2; d_xy; d_xz; d_y2; d_yz; d_z2];
                X=5;
            end     
            
            cart_exp = [cart_exp; subshells]; % mx3 array of the cartesian 
                                              % exponenets where m is the 
                                              % number of basis functions.

            for t=0:X
                rad_exp{q+t} = bdef{K}(i).exponents;
                coeffs{q+t} = bdef{K}(i).coeffs(1,:);
                if (t>=1) && (ifSP==2)
                    coeffs{q+t} = bdef{K}(i).coeffs(2,:);
                end
            end
            q=q+X+1;

        end
        
        [m,n] = size(cart_exp); % number of basis functions given by m.
                                % (number of rows in cartesian
                                % exponent array)
                                                               
        for p=1:m
            basis(p).atom = atoms(K);   % all m basis functions for atom K must be atom K
            basis(p).A = xyz_a0(K,:);   % same principle, but for nuclear coordinates
            basis(p).a = cart_exp(p,:); % each row in cart_exp corresponds to cartesian 
                                        % coordinate of the mth basis function
            
            basis(p).alpha = rad_exp{p};
            basis(p).d = coeffs(p);
            basis(p).N = (2/pi)^(3/4)*(2^sum(basis(p).a))*basis(p).alpha.^(((2*sum(basis(p).a))+3)/4)/sqrt(fac2(2*basis(p).a(1)-1)*fac2(2*basis(p).a(2)-1)*fac2(2*basis(p).a(3)-1));                          
        end      
          
    end
    
    %return basis;
    
end