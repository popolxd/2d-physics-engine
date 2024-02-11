# 2d physics engine

Collision are detected using the point line intersection. However, when one object is stationary, my implementation of the code doesn't work. Additionally sometimes it seems to just skip the collision even when both objects are moving.

Collision are resolved using impulse calculation, you can choose between elastic and inelastic collision. Only elastic ones work decently well for now. Most of the collision bugs happens here, but it doesn't mean that it is shit, I think that my implementation is correct.

todo:
- make the collisions viable even when one object is stationary
- make the collisions between bodies and other bodies (for now only player-body collision)
- try sat implementation (it is not that accurate, but it doesn't rely on floating point accuracy much)
- improve body-wall collisions
- make cross and dot product into macro functions
- optimize code