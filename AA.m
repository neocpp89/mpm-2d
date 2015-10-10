function []  = AA(fname, pstring)
    m = 0.219;
    v = 2:0.25:5.25;
    % v = 2:0.25:3;
    disp(v)
    ke = 0.5*m*v.^2;

    A = dlmread(fname);
    ke = ke(1:(size(A, 2)/2))
    lastpos = A(end, 2:2:end);
    firstpos = A(1, 2:2:end);
    depth = firstpos - lastpos;
    disp(depth)
    plot(ke, depth, pstring);
end
