# N-body Project

The only python file you will need to run for this project is nbody_fast.py. It can run any art with any boundary conditions
depending on th command-line arguments. You will need to specify the number of the part you want to run as thefirst argument,
then you can specify the boundary conditions to use using -BC p for persiodic or -BC np for non-periodic. By default the code
will run with periodic boundary condition, so you will only really need to specify them in part 3, for the non-periodic case.
The aimations are saved as mp4, so I don't think github will be ale to show those.

## Part 1

As can be seen in the part_1_p.mp4 video, a single particle stays at rest (the color of the dot changes at each iteration to 
show the this in not just a still image). I made sure that a particle's potential was symmetric from its positio to avoid 
a particle feeling its own potential. Havein periodic boundary conditions does not really influence this, because the potential
dies off well before reaching the particle by looping around the universe

## Part 2

As can be seen in part_2_p.mp4, two particles start in an approximately circular orbit, and stay in that orbit for a few 
turns. I fact, they stay in orbit for longer than what is shown in the video, but the orbits starts to become less circular.

## Part 3

For periodic boundary conditions, it seems the particles start to form what looks like spiral galaxies, then multiple of these 
galaxies merge into bigger, more slliptical clusters. This progression can be seen in part_3_p.mp4

For non-periodic boundary conditions, the process is initially very similar, but the particles all move towards the middle, and
instead of forming mutiple bigger galaxies, they all merge into one big cluster which pretty much instantly explodes and scatters
the particles all over our universe. This one can be seen in part_3_np.mp4

In both of these parts, energy changes quite a bit, meaning that there is room for improvement in the way I eveolve. Mostly,
I think the density grid could be found much more accurately than with the lazy method I used, where I just find the closest
grid cell to each mass. Interestingly, there is a spike in energy when the big cluster is formed in non-periodic BC.

## Part 4

Here, some structure appears, but a lot of the masses are so small that they seem to not really feel the potentialm so they 
essetially remain immobile during the entire simulation. Nonetheless, some clusters form, although they eventually breeak
down and go back to noise. Here, energy is pretty well conserved.
