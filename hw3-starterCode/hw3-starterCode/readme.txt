Assignment #3: Ray tracing

FULL NAME: Rachel Wang


MANDATORY FEATURES
------------------

<Under "Status" please indicate whether it has been implemented and is
functioning correctly.  If not, please explain the current status.>

Feature:                                 Status: finish? (yes/no)
-------------------------------------    -------------------------
1) Ray tracing triangles                  yes

2) Ray tracing sphere                     yes

3) Triangle Phong Shading                 yes

4) Sphere Phong Shading                   yes

5) Shadows rays                           yes

6) Still images                           yes
   
7) Extra Credit (up to 20 points)

   a) i implemented soft shadows by modeling each light source in a hemisphere,
   and choosing M random points on that hemisphere. i took the proportion of the
   light contribution from each of the "sub" lights, and added those to the final
   color of the pixel. (this creates M * num_lights amount of lights)

   b) i implemented antialiasing techniques by taking into account four rays per pixel
   (which i got from a stackexchange post), and averaging the pixel color of the four rays
   in my draw_scene() method. this works by supersampling the pixels.
