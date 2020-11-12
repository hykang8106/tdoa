function [] = learn_rng(use_rng)

if use_rng
    rng('default');
end

rand(1)

randi(10, 1)

rand(1)

randi(10, 1)

end