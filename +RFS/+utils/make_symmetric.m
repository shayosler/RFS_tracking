function out = make_symmetric(in)
%make_symmetric(in) Make a matrix that may be slightly asymmetric d/t floating
% point errors symmetric

out = (in + in') ./ 2;
end
