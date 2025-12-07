# RIO2 AVL Notes 

version 1

AVL 3.36

## Overview

From messing around with the current configuration in AVL, the plane overall just performs werll. I did some messing around with it like setting the plane to trim at various CL, which it didnt need much elevator deflection to do so (good sign). I did have issues rnning the rudder, but that's just my file being weird and isn't as critical as elevators and aileron control. I outputted the stability derivatives in txt files so I can look at them later, but the plane is stable and I just hope the yaw moment doesn't mess with the other derivatives. 

I did 3 different configurations, one for mission two and two for mission three. I'm not sure how we can simulate the shift in cg with avl, but with the derivatives there's probably something we can do with a state space simulation. 

## How to run

Compiling AVl may or may not be pain depending on your patience, but I'm running on 3.36 so if you decide to mess with this, I'd recomend the same AVL version. Documentation is very old school, but you can generally follow the steps below. I can send my workspace (minus the exe). 

```bash
cd AVL	# workspace
./avl	# run exe
load runs/rio2_m2.avl	# or whatever other case you want 
oper
g	# make sure the graphics work
d1 rm 0	# aileron to set roll to 0
d2 pm 0 # elevator to set pitch to 0
x	# run
t	# opern terretz plot
st	# hit enter after, shows the stability derivatives
``` 

## Next Steps 

I think we can close this for now, so controls are basically done yay. It may be worth getting a feel for AVL, but its not critical for our subteam. 

Next is most likely to focus on propeller selection, so maybe we can figure out how we can find our mission2 propeller and battery. Once we know the mission three propeller, we can use a 6S battery and get an optimal (and available) motor kv. Since both missions need the same motor, we can use that motor kv and simulate different propellers with different battery, between 3 and 12 S, and just go from there. 

