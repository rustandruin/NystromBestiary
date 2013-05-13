% out = dummy_extension(in)
%
% Returns all random numbers, for testing

function out = dummy_extension(in)

out.specerr = rand(2,in.q);
out.froerr = rand(2,in.q);
out.trerr = rand(2,in.q);
out.timings = rand(2,in.q);
out.numiters = rand(1,in.q);

end
