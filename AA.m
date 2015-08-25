function []  = AA(fname, pstring)
    m = 0.219;
    v = 1.5:0.25:3.75;
    ke = 0.5*m*v.^2;

    A = load(fname);
    lastpos = A(end, 2:2:end);
    firstpos = A(1, 2:2:end);
    depth = firstpos - lastpos;
    plot(ke, depth, pstring);
end
