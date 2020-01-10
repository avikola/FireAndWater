# Fire & Smoke Simulation

[Jump to Results](#results)

An Interactive OpenGL Smoke & Fire Simulation

Debugging: (Release x86)


## Fire

* Attempted basic rendering of fire with GL_POINTS for particles.
    * Set up decreasing life span with linearly changing colours for the particles
* Issue: Doesn't look realistic due to lack of texturing.

* Used fire texture on quads instead.
    * Implemented a shrinking effect for the diamond quads with decreasing life.


## Smoke

* Followed Jos Stam's work in Stable Fluids.
    * i.e. Navier Stoke's method with Mass Conservation Condition.
    
* 4 major components :
    * Forces - asserting forces on generated smoke in a given direction.
    * Advection - propogation of the fluid.
    * Diffusion - handling viscosity - our smoke didn't need to be viscous.
    * Projection - to maintain the mass conservation condition.
    

## UI

Various options to for users to interact with the simulation and control parameters:

* Mouse Interaction
    * Right-click hold to generate additional smoke at mouse location.
    * Left-click drag to distort existing smoke.

* Keyboard Interacion
    * Change position of smoke origination.
    * Change direction of smoke propogation.
    * Change velocity of smoke generation.
    * Adjust smoke density.
    * Change view mode between smoke, fire, fire + smoke.



## Results

### Smoke + UI


![](RealisticFire/RealisticFire/images/smoke_render.png)


### Basic GL_POINTS Fire


![](RealisticFire/RealisticFire/images/points_fire.png)


### Textured Fire


![](RealisticFire/RealisticFire/images/textured_fire.png)



### Fire + Smoke


![](RealisticFire/RealisticFire/images/fire_and_smoke.png)

<br/>
<br/>
<br/>
<br/>

<p align="center"><i>In Collaboration with <a href="https://github.com/mrobinson242">Matthew Robinson<a> and <a href="https://github.com/14joshb">Aditya Josh Batra</a></i></p>
