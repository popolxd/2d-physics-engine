# 2d physics engine

Collision are detected using the point line intersection. Recently I have fixed the bugs with collision detection. I seems to be correct now.

Collision are resolved using impulse calculation, you can choose between elastic and inelastic collision. But inelastic collisions are really unrealistic. Resolving of the collision is in my opinion correct, but it fails when objects are rotation at really high angular velocities, but I think it is becouse when I check the movements of the points, I assume it is linear, when in reality it is curved slightly. Higher angular velocity, more curvature.

todo:
- make the collisions between bodies and other bodies (for now only player-body collision)
- try to make point-line collision detection that also supports rotation
- make it possible for 2 collisions to occur in one frame
- make cross and dot product into macro functions
- optimize code