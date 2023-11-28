function h = pos_L1_hess(R, Z, nuk, U, data, lbd)

% compute RRZ
if isfield(data, 'RRZ')
    RRZ = data.RRZ;
else
    RRZ = R' * R - Z / nuk;
end

pgU = U'*R;
pgU = pgU + pgU';
pgU = pgU .*(abs(RRZ)<lbd/nuk);

% compute H2
h = R * pgU;
end