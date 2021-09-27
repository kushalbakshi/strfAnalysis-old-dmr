function [e_spectral, e_temporal, i_spectral, i_temporal] = findSTRFbw(STA, taxis, faxis)

excitatory = STA>0;
[r, c] = find(excitatory == 1);
c = unique(c);
c_dif = find(diff([1;c;1])~=1);
[lmax, jmax] = max(diff(c_dif));
start_max = c(jmax);
e_temporal = range(taxis(start_max:lmax));

r = unique(r);
r_dif = find(diff([1;r;1])~=1);
[lmax, jmax] = max(diff(r_dif));
end_max = r(jmax);
e_spectral = range(faxis(1, lmax:end_max)/1000);

inhibitory = STA < 0;
[r, c] = find(inhibitory == 1);
c = unique(c);
c_dif = find(diff([1;c;1])~=1);
[lmax, jmax] = max(diff(c_dif));
start_max = c(jmax);
i_temporal = range(taxis(start_max:lmax));

r = unique(r);
r_dif = find(diff([1;r;1])~=1);
[lmax, jmax] = max(diff(r_dif));
end_max = r(jmax);
i_spectral = range(faxis(1, lmax:end_max)/1000);