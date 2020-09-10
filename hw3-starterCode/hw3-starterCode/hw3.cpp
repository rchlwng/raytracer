/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Rachel Wang
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

#define WIDTH 640
#define HEIGHT 480
float aspectRatio = (float)WIDTH / (float)HEIGHT;
float zValue = -1;
float smallConstant = .00001;
float lightRadius = .05;


//the field of view of the camera
#define fov 60.0
#define M 50

#define PI 3.14159265

unsigned char buffer[HEIGHT][WIDTH][3];

struct Intersection {
  bool intersects;
  glm::vec3 intersectionPoint;
};

struct Ray {
  glm::vec3 origin;
  glm::vec3 direction;
};

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

float screenScale = 1.0;
float tangentMultiplier = (float)tan((fov / 2) * PI/180.0f);

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

glm::vec3 clamp(glm::vec3 rawValues) {
  if (rawValues.x > 1.0) {
    rawValues.x = 1.0;
  }
  else if (rawValues.x < 0) {
    rawValues.x = 0;
  }

  if (rawValues.y > 1.0) {
    rawValues.y = 1.0;
  }
  else if (rawValues.y < 0) {
    rawValues.y = 0;
  }

  if (rawValues.z > 1.0) {
    rawValues.z = 1.0;
  }
  else if (rawValues.z < 0) {
    rawValues.z = 0;
  }
  return rawValues;

} 

// stackoverflow
double randValue(double min, double max) {
  double f = (double)rand() / RAND_MAX;
  return min + f * (max - min);
}

// stackoverflow, marsaglia's algorithm
std::vector<Light> circularLight(Light light) {

  std::vector<Light> lightSamples;

  // stackeroverflow, randomly sampling a hemisphere
  for (int i = 0; i < M; i++) {

    double x1 = randValue(-1, 1);
    double x2 = randValue(-1, 1);

    if (pow(x1, 2) + pow(x2, 2) >= 1) {
      x1 = 0;
      x2 = 0;
    }

    // spherical coordinates diagram (cartesian conversion)
    double x = 2 * x1 * sqrt(1 - pow(x1, 2) - pow(x2, 2)) + light.position[0];
    double y = 2 * x2 * sqrt(1 - pow(x1, 2) - pow(x2, 2)) + light.position[1];
    double z = 1 - 2 * (pow(x1, 2) + pow(x2, 2)) + light.position[2];

    Light sampledLight;
    sampledLight.position[0] = x;
    sampledLight.position[1] = y;
    sampledLight.position[2] = z;
    sampledLight.color[0] = (double)light.color[0]/(double)M;
    sampledLight.color[1] = (double)light.color[1]/(double)M;
    sampledLight.color[2] = (double)light.color[2]/(double)M;

    lightSamples.push_back(sampledLight);

  }

  return lightSamples;

}

std::vector<Ray> shootRay(float x, float y) {
  // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-generating-camera-rays/generating-camera-rays

  // normalized device coordinates
  float preliminaryXOne = (x + .25) / WIDTH;
  float preliminaryYOne = (y + .25) / HEIGHT;

  // screen coords
  float screenXOne = 2.0 * preliminaryXOne - 1;
  float screenYOne = 2.0 * preliminaryYOne - 1; /*2.0 * preliminaryY - 1;*/

  float xDestinationOne = screenXOne * aspectRatio * tangentMultiplier;
  float yDestinationOne = screenYOne * tangentMultiplier;

  // powerpoint 15, slide 10
  glm::vec3 newRayOne = glm::normalize(glm::vec3(xDestinationOne, yDestinationOne, zValue));
  Ray rayOne = { glm::vec3(0, 0, 0), newRayOne };

  // normalized device coordinates
  float preliminaryXTwo = (x + .25) / WIDTH;
  float preliminaryYTwo = (y + .75) / HEIGHT;

  // screen coords
  float screenXTwo = 2.0 * preliminaryXTwo - 1;
  float screenYTwo = 2.0 * preliminaryYTwo - 1; /*2.0 * preliminaryY - 1;*/

  float xDestinationTwo = screenXTwo * aspectRatio * tangentMultiplier;
  float yDestinationTwo = screenYTwo * tangentMultiplier;

  // powerpoint 15, slide 10
  glm::vec3 newRayTwo = glm::normalize(glm::vec3(xDestinationTwo, yDestinationTwo, zValue));
  Ray rayTwo = { glm::vec3(0, 0, 0), newRayTwo };

  // normalized device coordinates
  float preliminaryXThree = (x + .75) / WIDTH;
  float preliminaryYThree = (y + .75) / HEIGHT;

  // screen coords
  float screenXThree = 2.0 * preliminaryXThree - 1;
  float screenYThree = 2.0 * preliminaryYThree - 1; /*2.0 * preliminaryY - 1;*/

  float xDestinationThree = screenXThree * aspectRatio * tangentMultiplier;
  float yDestinationThree = screenYThree * tangentMultiplier;

  // powerpoint 15, slide 10
  glm::vec3 newRayThree = glm::normalize(glm::vec3(xDestinationThree, yDestinationThree, zValue));
  Ray rayThree = { glm::vec3(0, 0, 0), newRayThree };

  // normalized device coordinates
  float preliminaryXFour = (x + .75) / WIDTH;
  float preliminaryYFour = (y + .25) / HEIGHT;

  // screen coords
  float screenXFour = 2.0 * preliminaryXFour - 1;
  float screenYFour = 2.0 * preliminaryYFour - 1; /*2.0 * preliminaryY - 1;*/

  float xDestinationFour = screenXFour * aspectRatio * tangentMultiplier;
  float yDestinationFour = screenYFour * tangentMultiplier;

  // powerpoint 15, slide 10
  glm::vec3 newRayFour = glm::normalize(glm::vec3(xDestinationFour, yDestinationFour, zValue));
  Ray rayFour = { glm::vec3(0, 0, 0), newRayFour };

  std::vector<Ray> antialiased;
  antialiased.push_back(rayOne);
  antialiased.push_back(rayTwo);
  antialiased.push_back(rayThree);
  antialiased.push_back(rayFour);

  return antialiased;

}

glm::vec3 reflection(glm::vec3 normal, glm::vec3 light) {

  float ln = glm::dot(light, normal);
  glm::vec3 reflectedRay = (2 * ln * normal) - light;
  reflectedRay = glm::normalize(reflectedRay);
  return reflectedRay;

}

glm::vec3 phongSphere(glm::vec3 intersection, Sphere sphere, Light light) {

  // phong equation: L(kd(l*n) + ks(r*v)^a)
  glm::vec3 normal = intersection - glm::vec3(sphere.position[0],
   sphere.position[1], sphere.position[2]);

  normal *= 1 / sphere.radius;
  glm::vec3 lightDir = glm::vec3(light.position[0], light.position[1], light.position[2])
   - intersection;
  lightDir = glm::normalize(lightDir);
  glm::vec3 reflectedRay = reflection(normal, lightDir);
  glm::vec3 viewer = intersection * -1.0f;
  viewer = glm::normalize(viewer);

  float ln = glm::dot(lightDir, normal);
  float rv = glm::dot(reflectedRay, viewer);

  if (ln > 1.0) {
    ln = 1.0;
  }
  else if (ln < 0) {
    ln = 0;
  }

  if (rv > 1.0) {
    rv = 1.0;
  }
  else if (rv < 0) {
    rv = 0;
  }

  glm::vec3 pixelColor;

  // rgb
  pixelColor.x = light.color[0] * (sphere.color_diffuse[0] * ln + (sphere.color_specular[0] * pow(rv, sphere.shininess)));
  pixelColor.y = light.color[1] * (sphere.color_diffuse[1] * ln + (sphere.color_specular[1] * pow(rv, sphere.shininess)));
  pixelColor.z = light.color[2] * (sphere.color_diffuse[2] * ln + (sphere.color_specular[2] * pow(rv, sphere.shininess)));

  pixelColor = clamp(pixelColor);
  return pixelColor;

}

glm::vec3 phongTriangle(glm::vec3 intersection, Triangle triangle, Light light) {

  // powerpoint 16, slide 10
  glm::vec3 a = glm::vec3(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
  glm::vec3 b = glm::vec3(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
  glm::vec3 c = glm::vec3(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);

  // thank you, stackexchange post :)
  glm::vec3 sideOne = b - a;
  glm::vec3 sideTwo = c - b;
  glm::vec3 sideThree = a - c;

  glm::vec3 normalEdge = c - a;
  glm::vec3 planeNormal = glm::cross(sideOne, normalEdge);

  // stackexchange
  float areaWhole = glm::dot(planeNormal, glm::cross(sideOne, normalEdge));
  float areaBCIntersect = glm::dot(planeNormal, glm::cross(b-intersection, c-intersection));
  float areaCAIntersect = glm::dot(planeNormal, glm::cross(c-intersection, a-intersection));

  float alpha = areaBCIntersect / areaWhole;
  float beta = areaCAIntersect / areaWhole;
  float gamma = 1 - alpha - beta;

  // origin point a

  float normalX = alpha * triangle.v[0].normal[0] + beta * triangle.v[1].normal[0] + gamma * triangle.v[2].normal[0];
  float normalY = alpha * triangle.v[0].normal[1] + beta * triangle.v[1].normal[1] + gamma * triangle.v[2].normal[1];
  float normalZ = alpha * triangle.v[0].normal[2] + beta * triangle.v[1].normal[2] + gamma * triangle.v[2].normal[2];

  glm::vec3 normal = glm::vec3(normalX, normalY, normalZ);
  normal = glm::normalize(normal);
  glm::vec3 lightDir = glm::vec3(light.position[0], light.position[1], light.position[2]) - intersection;
  lightDir = glm::normalize(lightDir);
  glm::vec3 reflectedRay = reflection(normal, lightDir);
  glm::vec3 viewer = intersection * (float)-1;
  viewer = glm::normalize(viewer);

  float ln = glm::dot(lightDir, normal);
  float rv = glm::dot(reflectedRay, viewer);

  if (ln > 1.0) {
    ln = 1.0;
  }
  else if (ln < 0) {
    ln = 0;
  }

  if (rv > 1.0) {
    rv = 1.0;
  }
  else if (rv < 0) {
    rv = 0;
  }


  // interpolating shader values
  float diffuseR = alpha * triangle.v[0].color_diffuse[0] + beta * triangle.v[1].color_diffuse[0] + gamma * triangle.v[2].color_diffuse[0];
  float diffuseG = alpha * triangle.v[0].color_diffuse[1] + beta * triangle.v[1].color_diffuse[1] + gamma * triangle.v[2].color_diffuse[1];
  float diffuseB = alpha * triangle.v[0].color_diffuse[2] + beta * triangle.v[1].color_diffuse[2] + gamma * triangle.v[2].color_diffuse[2];

  float specularR = alpha * triangle.v[0].color_specular[0] + beta * triangle.v[1].color_specular[0] + gamma * triangle.v[2].color_specular[0];
  float specularG = alpha * triangle.v[0].color_specular[1] + beta * triangle.v[1].color_specular[1] + gamma * triangle.v[2].color_specular[1];
  float specularB = alpha * triangle.v[0].color_specular[2] + beta * triangle.v[1].color_specular[2] + gamma * triangle.v[2].color_specular[2];

  float shininess = alpha * triangle.v[0].shininess + beta * triangle.v[1].shininess + gamma * triangle.v[2].shininess;

  glm::vec3 pixelColor;

  // rgb
  pixelColor.x = light.color[0] * (diffuseR * ln + specularR * pow(rv, shininess));
  pixelColor.y = light.color[1] * (diffuseG * ln + specularG * pow(rv, shininess));
  pixelColor.z = light.color[2] * (diffuseB * ln + specularB * pow(rv, shininess));

  pixelColor = clamp(pixelColor);

  return pixelColor;

}

Intersection sphereIntersect(Ray ray, Sphere sphere) {
  
  // powerpoint 16, slide 6
  float centerX = sphere.position[0];
  float centerY = sphere.position[1];
  float centerZ = sphere.position[2];

  glm::vec3 spherePosition(centerX, centerY, centerZ);
  glm::vec3 origin = ray.origin - spherePosition;

  float b = 2.0 * (glm::dot(ray.direction, origin));
  float c = glm::dot(origin, origin) - pow(sphere.radius, 2);

  float discriminant = pow(b, 2) - (4 * c);

  Intersection none = { false, glm::vec3(0, 0, 0) };

  float t0, t1;

  // no solutions to quadratic equation
  if (discriminant < 0) {
    return none;
  }
  // one real root
  else if (discriminant == 0) {
    t0 = -b / 2.0f;
    t1 = -b / 2.0f;
  }
  else {
    t0 = (-b + sqrt(discriminant)) / 2.0;
    t1 = (-b - sqrt(discriminant)) / 2.0;
  }

  if (t0 < smallConstant || t1 < smallConstant) {
    return none;
  }
  else {

    if (t0 > t1) {
      std::swap(t0, t1);
    }

    float finalT = t0;
    glm::vec3 intersectionRay = ray.origin + (ray.direction * finalT);
    Intersection success = { true, intersectionRay };
    return success;
  }

}

Intersection triangleIntersect(Ray ray, Triangle triangle) {

  // powerpoint 16, slide 10
  glm::vec3 a = glm::vec3(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
  glm::vec3 b = glm::vec3(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
  glm::vec3 c = glm::vec3(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);

  // thank you, stackexchange post :)
  glm::vec3 sideOne = b - a;
  glm::vec3 sideTwo = c - b;
  glm::vec3 sideThree = a - c;

  glm::vec3 normalEdge = c - a;
  // origin point a
  glm::vec3 normal = glm::cross(sideOne, normalEdge);
  normal = glm::normalize(normal);


  // powerpoint 16, slide 11
  float nD = glm::dot(normal, ray.direction);
  // d constant in plane equation
  float d = glm::dot(normal, ray.origin-a);
  float t = -1;
  Intersection none = { false, glm::vec3(0, 0, 0) };

  if (nD != 0) {
    t = -d / nD;
  }

  if (nD == 0 || t <= smallConstant) {
    return none;
  }
  else {
    
    // barycentric coordinates
    // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates
    
    glm::vec3 possibleIntersection = ray.origin + (ray.direction * t);
    glm::vec3 orthogonalVecA = glm::cross(sideOne, possibleIntersection - a);
    glm::vec3 orthogonalVecB = glm::cross(sideTwo, possibleIntersection - b);
    glm::vec3 orthogonalVecC = glm::cross(sideThree, possibleIntersection - c);
    
    // ray does not hit the triangle
    if (glm::dot(normal, orthogonalVecA) < 0 || 
     glm::dot(normal, orthogonalVecB) < 0 ||
      glm::dot(normal, orthogonalVecC) < 0) {
      return none;
    }
    else {
      Intersection success = { true, ray.origin + (ray.direction * t) };
      return success;
    }
    
  }

}


glm::vec3 rayTrace(Ray ray) {

  // default to white
  glm::vec3 color = glm::vec3(.24, .23, .5);

  float closestSphere = -pow(10, 10);
  glm::vec3 closestSpherePoint;
  int closestSphereIndex = -1;

  for (int i = 0; i < num_spheres; i++) {
    Intersection intersection = sphereIntersect(ray, spheres[i]);
    if (intersection.intersects &&
     intersection.intersectionPoint.z > closestSphere) {
       closestSphere = intersection.intersectionPoint.z;
       closestSpherePoint = intersection.intersectionPoint;
       closestSphereIndex = i;
    }
  }

  float closestTriangle = -pow(10, 10);
  glm::vec3 closestTrianglePoint;
  int closestTriangleIndex = -1;
  for (int i = 0; i < num_triangles; i++) {
    Intersection intersection = triangleIntersect(ray, triangles[i]);
    if (intersection.intersects &&
     intersection.intersectionPoint.z > closestTriangle) {
       closestTriangle = intersection.intersectionPoint.z;
       closestTrianglePoint = intersection.intersectionPoint;
       closestTriangleIndex = i;
    }
  }

  // powerpoint 15, slide 22
  if (closestTriangleIndex != -1 || closestSphereIndex != -1) {
    color = glm::vec3(0, 0, 0);
    // sphere is closer
    if (closestSphere > closestTriangle) {

      // for each light source
      for (int i = 0; i < num_lights; i++) {

        std::vector<Light> lightSamples = circularLight(lights[i]);

        for (int z = 0; z < M; z++) {

            bool blocked = false;
            glm::vec3 lightSource = glm::vec3(lightSamples[z].position[0],
            lightSamples[z].position[1], lightSamples[z].position[2]);
            Ray shadowRay = { closestSpherePoint, glm::normalize(lightSource - closestSpherePoint) };

            // does it hit an opaque object (spheres)
            for (int k = 0; k < num_spheres; k++) {

                Intersection intersectionSphere = sphereIntersect(shadowRay, spheres[k]);

                // not the current sphere
                if (intersectionSphere.intersects && k != closestSphereIndex) {

                  float distanceFromLight = glm::distance(lightSource, closestSpherePoint);
                  float distanceFromSphere = glm::distance(intersectionSphere.intersectionPoint, closestSpherePoint);

                  // light source is blocked by errant sphere
                  if (distanceFromSphere < distanceFromLight) {
                    blocked = true;
                    break; // no need to keep iterating
                  }

                }
            }
            // does it hit an opaque object (triangles)
            // FIXED BUG: comparing distances of wrong objects
            for (int j = 0; j < num_triangles; j++) {

              if (blocked) {
                break;
              }

              Intersection intersectionTri = triangleIntersect(shadowRay, triangles[j]);
              if (intersectionTri.intersects) {

                float distanceFromLight = glm::distance(lightSource, closestSpherePoint);
                float distanceFromTri = glm::distance(intersectionTri.intersectionPoint, closestSpherePoint);

                if (distanceFromTri < distanceFromLight) {
                  blocked = true;
                  break; // no need to keep iterating
                }
              }
            }
            if (!blocked) {
              glm::vec3 temp = clamp(phongSphere(closestSpherePoint, spheres[closestSphereIndex], lightSamples[z]));
              color += temp;
            }

        }
        
      }
    }
    // triangle
    else if (closestTriangle > closestSphere) {
      for (int i = 0; i < num_lights; i++) {

        std::vector<Light> lightSamples = circularLight(lights[i]);
        for (int z = 0; z < M; z++) {

          bool blocked = false;
          glm::vec3 lightSource = glm::vec3(lightSamples[z].position[0],
            lightSamples[z].position[1], lightSamples[z].position[2]);
          Ray shadowRay = { closestTrianglePoint, glm::normalize(lightSource - closestTrianglePoint) };

          // checking for opaque objects again
          for (int k = 0; k < num_spheres; k++) {

              Intersection intersectionSphere = sphereIntersect(shadowRay, spheres[k]);
              if (intersectionSphere.intersects) {

                float distanceFromLight = glm::distance(lightSource, closestTrianglePoint);
                float distanceFromSphere = glm::distance(intersectionSphere.intersectionPoint, closestTrianglePoint);

                if (distanceFromSphere < distanceFromLight) {
                  blocked = true;
                  break; // no need to keep iterating
                }

              }
          }
          for (int j = 0; j < num_triangles; j++) {

            if (blocked) {
              break;
            }

            Intersection intersectionTri = triangleIntersect(shadowRay, triangles[j]);
            if (intersectionTri.intersects && j != closestTriangleIndex) {

              float distanceFromLight = glm::distance(lightSource, closestTrianglePoint);
              float distanceFromTri = glm::distance(intersectionTri.intersectionPoint, closestTrianglePoint);

              if (distanceFromTri < distanceFromLight) {
                blocked = true;
                break; // no need to keep iterating
              }
            }
          }
          
          if (!blocked) {
            glm::vec3 temp = clamp(phongTriangle(closestTrianglePoint, triangles[closestTriangleIndex], lightSamples[z]));
            color += temp;
          }

        }
      }
    }
  }

  glm::vec3 ambientLight = glm::vec3(ambient_light[0], ambient_light[1], ambient_light[2]);
  color += ambientLight;

  color = clamp(color);

  return color;

}

//MODIFY THIS FUNCTION
void draw_scene()
{

  // backwards ray tracing
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      // stackexchange antialiasing, four subpixels per pixel
      std::vector<Ray> antialiased = shootRay(x, y);
      glm::vec3 color = glm::vec3(0, 0, 0);
      for (int i = 0; i < 4; i++) {
        glm::vec3 pixelColor = rayTrace(antialiased[i]);
        color += pixelColor;
      }
      color /= (float)4;
      color = clamp(color);
      plot_pixel(x, y, color.x * 255, color.y * 255, color.z * 255);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

