# Cislunar-Force-Model
This assembly of codes in Matlab allows the propagation of a satellite in the solar system with very high accuracy. The emphasis is made on cislunar and circumlunar spacecraft cases.  

1) First of all, this Matlab Code call JPL SPICE ephemeris data for the position/velocity of other celestial bodies in order to have the ephemeris precision level. Therefore you will need to have downloaded MICE (the Matlab SPICE) to use this code. You will also have to download the ephemeris for the Moon etc You can find everything here: https://naif.jpl.nasa.gov/naif/toolkit.html
2) This code is a high-fidelity force model made specifically for cislunar and circumlunar spacecraft. However, it is totally possible to modify it for other applications around other central bodies.
3) The following forces can be taken into account depending on the user's propagation objectives:
  - 2BP + Non-sphericity of the central body using Gravity Spherical Harmonics (GSH) under the Gottlieb formulation (to avoid any singularity and keep a good precision with high order and degree). Some documents are attached in the main branch of this project in order to give the interested reader further details about the derivation of such formulation as well as its implementation which require many steps and is mathematically involved
  - 3rd body. Any third celestial body from the solar system (thanks to SPICE ephemeris datasets). However, for some of the,m you might need to download the specific ephemeris that cover the requested 3rd body during the required time period of your propagation
  - Solar Radiation Pressure (SRP) from the photon of the Sun. The spacecraft (S/C) is modeled with the sun's rays hitting the S/C perpendicularly. The area of the surface facing the Sun will be needed. The user has to choose the shadow model: conical (higher fidelity), cylindrical (lower fidelity), or none (lowest fidelity).
  - Solid Tides. A bit more advanced effects in astrodynamics formulations but quite significant for the Earth-Moon system. If studying a central body with water ocean, might want to add the Ocean tides effects. In any case, these effects are taken into account via a modification of the GSH coefficients (see the references for more details)
  - Relativity effects expressed in a Newtonian formulation. Multiple effects can be taken into account (Schwarzschild correction and Lense-Thirring correction)
4) Keep the length of the propagation quite small as it is using a lot of computer resources due to the ephemeris precision data
5) The outcome of the propagated trajectories has been compared to GMAT and gives very good results. (less than 1m over 15 days of propagation for low-altitude orbiters that are therefore in a highly perturbed environment).

If you have any questions/remarks, please contact me, I would be happy to answer you!
