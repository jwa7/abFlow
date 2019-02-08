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

function buildbasis(atoms,xyz_a0,basissetdef)
    s =    [0 0 0];
    p_x =  [1 0 0];
    p_y =  [0 1 0];
    p_z =  [0 0 1];
    d_x2 = [2 0 0];
    d_xy = [1 1 0];
    d_xz = [1 0 1];
    d_y2 = [0 2 0];
    d_yz = [0 1 1];
    d_z2 = [0 0 2];
    
    basis = struct;
    
    for p=1:numel(atoms) % iterating over each atom in the atom list of the molecule
       
        basis(p).atom = atoms(p);
        basis(p).A = xyz_a0(p,:);
        
        basis(p).a = [];
        
        for i=1:numel(basissetdef{p}(i).shelltype) % iterating over 
            
            type = basissetdef{p}(i).shelltype;
            
            if type == 'S'
                cart_exp = [s];
            elseif type == 'SP'
                cart_exp = [s; p_x; p_y; p_z];
            elseif type == 'P'
                cart_exp = [p_x; p_y; p_z];
            elseif type == 'D'
                cart_exp = [d_x2; d_xy; d_xz; d_y2; d_yz; d_z2];
            end     
            
            basis(p).a = [basis(p).a; cart_exp];
            
        end
            
        
        shelltype = basissetdef{p}.shelltype;
        exponents = basissetdef{p}.exponents;
        coeffs = basissetdef{p}.coeffs;
        
    end



    return basis;
end