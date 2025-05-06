function [c] = MultiPermuations(m,r)
% a_1 + ... + a_r = m                                                                                                                                                                            
                                                                               
% All possible placements of internal dividers.                                       
dividers = nchoosek(1:(m+r-1), r-1);                                                  
ndividers = size(dividers, 1);                                                        
% Add dividers at the beginning and end.                                              
b = cat(2, zeros(ndividers, 1), dividers, (m+r)*ones(ndividers, 1));                  
% Find distances between dividers.                                                    
c = diff(b, 1, 2) - 1;

% num = factorial(m+r-1)/factorial(r-1)/factorial(m);
end