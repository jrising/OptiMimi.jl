func(x) = sin(2pi * x)

fs, mins, maxs = univariatematrix(func, -.1, .25, 0., .1)

linprog(fs, [], '<', [], mins, maxs)
