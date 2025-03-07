function coord = convert_coord(i,j,k,N1,N2,N3)
    coord = i + (N1+1)*(j-1) + (N1+1)*(N2+1)*(k-1);
end