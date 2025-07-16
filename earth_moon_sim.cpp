#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>

const double G   = 6.67430e-11;
const double M_s = 1.989e30;         // Sun mass (kg)
const double M_e = 5.972e24;         // Earth mass (kg)
const double M_m = 7.3477e22;        // Moon mass (kg)
const double r_se= 1.496e11;         // Sun–Earth distance (m)
const double r_em= 384400e3;         // Earth–Moon distance (m)
const double T_e = 365.25 * 24*3600; // Earth orbital period (s)
const double T_m = 27.322 * 24*3600; // Moon orbital period (s)
const double omega_e = 2*M_PI / T_e;
const double omega_m = 2*M_PI / T_m;

const double dt    = 60;              // simulation timestep (s)
const double t_max = 5 * 24*3600;     // run for 5 days

// Visualization scale (meters to pixels)
const double SCALE = 2.5e-9;

struct Vec2 {
    double x, y;
    Vec2 operator+(const Vec2& v) const { return {x+v.x, y+v.y}; }
    Vec2 operator-(const Vec2& v) const { return {x-v.x, y-v.y}; }
    Vec2 operator*(double s) const { return {x*s, y*s}; }
    Vec2 operator/(double s) const { return {x/s, y/s}; }
    double mag() const { return std::sqrt(x*x + y*y); }
};

Vec2 sun_pos() {
    return {0, 0};
}

Vec2 earth_pos(double t) {
    return {r_se * std::cos(omega_e * t), r_se * std::sin(omega_e * t)};
}

Vec2 moon_pos(double t) {
    Vec2 e = earth_pos(t);
    return e + Vec2{r_em * std::cos(omega_m * t), r_em * std::sin(omega_m * t)};
}

Vec2 acceleration(const Vec2& r, double t) {
    // Sun gravity
    Vec2 a_s = (sun_pos() - r) * (G * M_s / std::pow((sun_pos() - r).mag(), 3));
    // Earth gravity
    Vec2 e = earth_pos(t);
    Vec2 a_e = (e - r) * (G * M_e / std::pow((e - r).mag(), 3));
    // Moon gravity
    Vec2 m = moon_pos(t);
    Vec2 a_m = (m - r) * (G * M_m / std::pow((m - r).mag(), 3));
    return a_s + a_e + a_m;
}

int main() {
    sf::RenderWindow window(sf::VideoMode(1200, 900), "Sun-Earth-Moon-Spacecraft Simulation");
    window.setFramerateLimit(60);

    // Initial spacecraft state (200 km above Earth, in Earth's direction)
    double t = 0.0;
    Vec2 r = earth_pos(0) + Vec2{6671e3, 0};
    Vec2 v = {0, 7800};

    // For trails
    sf::VertexArray craftTrail(sf::LineStrip);

    while (window.isOpen() && t < t_max) {
        sf::Event event;
        while (window.pollEvent(event))
            if (event.type == sf::Event::Closed) window.close();

        // RK4 integration
        Vec2 k1v = acceleration(r, t) * dt;
        Vec2 k1r = v * dt;

        Vec2 k2v = acceleration(r + k1r * 0.5, t + 0.5*dt) * dt;
        Vec2 k2r = (v + k1v * 0.5) * dt;

        Vec2 k3v = acceleration(r + k2r * 0.5, t + 0.5*dt) * dt;
        Vec2 k3r = (v + k2v * 0.5) * dt;

        Vec2 k4v = acceleration(r + k3r, t + dt) * dt;
        Vec2 k4r = (v + k3v) * dt;

        v = v + (k1v + k2v*2 + k3v*2 + k4v) / 6.0;
        r = r + (k1r + k2r*2 + k3r*2 + k4r) / 6.0;
        t += dt;

        // Draw
        window.clear(sf::Color::Black);

        // Sun
        sf::CircleShape sun(20);
        sun.setFillColor(sf::Color::Yellow);
        sun.setOrigin(20, 20);
        sun.setPosition(600, 450);
        window.draw(sun);

        // Earth
        Vec2 e = earth_pos(t);
        sf::CircleShape earth(8);
        earth.setFillColor(sf::Color::Blue);
        earth.setOrigin(8, 8);
        earth.setPosition(600 + e.x * SCALE, 450 + e.y * SCALE);
        window.draw(earth);

        // Moon
        Vec2 m = moon_pos(t);
        sf::CircleShape moon(4);
        moon.setFillColor(sf::Color(200, 200, 200));
        moon.setOrigin(4, 4);
        moon.setPosition(600 + m.x * SCALE, 450 + m.y * SCALE);
        window.draw(moon);

        // Spacecraft
        sf::CircleShape craft(3);
        craft.setFillColor(sf::Color::Red);
        craft.setOrigin(3, 3);
        craft.setPosition(600 + r.x * SCALE, 450 + r.y * SCALE);
        window.draw(craft);

        // Trail
        craftTrail.append(sf::Vertex(sf::Vector2f(600 + r.x * SCALE, 450 + r.y * SCALE), sf::Color::Green));
        window.draw(craftTrail);

        window.display();
    }
    return 0;
}