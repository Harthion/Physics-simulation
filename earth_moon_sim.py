from vpython import sphere, vector, rate, color, scene, textures
import numpy as np

#constants
G   = 6.67430e-11          # gravitational constant
M_e = 5.972e24             # Earth mass (kg)
M_m = 7.3477e22            # Moon mass  (kg)
r_em= 384400e3             # Earth–Moon distance (m)
T_m = 27.322 * 24*3600     # Moon orbital period (s)
ω   = 2*np.pi / T_m        # Moon angular speed

#scene setup
init_pos = vector(6671e3, 0, 0)   # start 200 km above Earth
init_vel = vector(0, 7800, 0)     # approx LEO speed
dt       = 60                     # simulation timestep (s)
t_max    = 5 * 24*3600            # run for 5 days

scene.width  = 800
scene.height = 600
scene.title  = "Earth–Moon–Spacecraft Simulation"


#earth eith texture
earth = sphere(
    pos=vector(0,0,0),
    radius=6.371e6,
    texture=textures.earth,
    shininess=0.8
)

# Starfield background
star_layer = box(
    pos=vector(0, 0, -1e9),
    size=vector(3e9, 3e9, 1e6),
    texture="https://upload.wikimedia.org/wikipedia/commons/0/01/Starsinthesky.jpg",
    emissive=True,
    opacity=1
)

moon = sphere(
    pos=vector(r_em,0,0),
    radius=1.737e6,
    color=color.white,
    make_trail=True,
    trail_radius=1e6
)

#spacecreaft with trail 
craft = sphere(
    pos=init_pos,
    radius=2e6,
    color=color.red,
    make_trail=True,
    trail_radius=5e5
)


def moon_pos(t):
    return vector(r_em*np.cos(ω*t),
                  r_em*np.sin(ω*t),
                  0)

def acceleration(r, t):
    # earthgrav
    a_e = G * M_e * (-r) / (r.mag**3)
    # moongrav
    r_m = moon_pos(t) - r
    a_m = G * M_m * r_m / (r_m.mag**3)
    return a_e + a_m

# AI below - inital state 
r = init_pos
v = init_vel
t = 0.0

#simulation loop 
while t < t_max:
    rate(100)   # cap to 100 iterations/sec

    # move Moon
    moon.pos = moon_pos(t)

    # RK4 integration for the craft
    k1v = acceleration(r, t) * dt
    k1r = v * dt

    k2v = acceleration(r + 0.5*k1r, t + 0.5*dt) * dt
    k2r = (v + 0.5*k1v) * dt

    k3v = acceleration(r + 0.5*k2r, t + 0.5*dt) * dt
    k3r = (v + 0.5*k2v) * dt

    k4v = acceleration(r + k3r, t + dt) * dt
    k4r = (v + k3v) * dt

    v = v + (k1v + 2*k2v + 2*k3v + k4v) / 6
    r = r + (k1r + 2*k2r + 2*k3r + k4r) / 6

    craft.pos = r
    t += dt
