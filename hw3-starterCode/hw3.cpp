/* **************************
 * CSCI 420
 * Assignment 3: Raytracer
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
#include "Math.h"
#include <math.h>
#include <iostream>

#define MAX_TRIANGLES 1000000
#define MAX_SPHERES 100000
#define MAX_LIGHTS 10000

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

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

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

struct Ray
{
    Vector3 startPos;
    Vector3 dir;
};

bool IntersectedSphere(Ray ray, double& t, int& sphereIndex)
{
    bool isIntersected = false;
    // for every sphere
    for (int i = 0; i < num_spheres; i++) {
        Vector3 sphere(
            spheres[i].position[0],
            spheres[i].position[1],
            spheres[i].position[2]
        );
        double b = -2.0 * Vector3::Dot(ray.dir,sphere);
        double c = Vector3::Dot(sphere, sphere) - pow(spheres[i].radius,2.0);
        double bSqaureMinus4c = pow(b,2.0) - 4 * c;
        if (bSqaureMinus4c < 0.0) { // abort if negative
            continue;
        }

        double t0 = (-b + sqrt(bSqaureMinus4c)) / 2.0;
        double t1 = (-b - sqrt(bSqaureMinus4c)) / 2.0;
        double tMin;
        if (t0 < 0.0 && t1 < 0.0) {
            continue;
        }
        else if (t0 < 0.0 && t1 > 0.0) {
            tMin = t1;
        }
        else if (t0 > 0.0 && t1 < 0.0) {
            tMin = t0;
        }
        else {
            tMin = Math::Min(t0, t1);
        }
        isIntersected = true;
        // closer intersection
        if (tMin < t) {
            t = tMin;
            sphereIndex = i;
        }
    }
    return isIntersected;
}

Vector3 GetTriangleNormal(Vector3 a, Vector3 b, Vector3 c)
{
    Vector3 normal = Vector3::Cross(b - a, c - a);
    return normal;
}

bool isPointInTri(Ray ray, double t, double& alpha, double& beta, double& gamma, Vector3& a, Vector3& b, Vector3& c)
{
    Vector3 triNormal = GetTriangleNormal(a, b, c);
    double tempA, tempB, tempG;
    // intersection point
    Vector3 p = ray.startPos + (ray.dir * t);

    // edges
    Vector3 ab = b - a;
    Vector3 bc = c - b;
    Vector3 ca = a - c;

    Vector3 ap = p - a;
    Vector3 bp = p - b;
    Vector3 cp = p - c;

    // whole triangle
    double whole = Vector3::Dot(triNormal, triNormal);
    if (whole == 0.0) // division by zero
    {
        return false;
    }

    Vector3 abXap = Vector3::Cross(ab, ap);
    if (Vector3::Dot(triNormal, abXap) < 0.0) {
        return false;
    }


    Vector3 bcXbp = Vector3::Cross(bc, bp);
    tempA = Vector3::Dot(triNormal, bcXbp) / whole;
    if (tempA < 0.0) {
        return false;
    }
    
    Vector3 caXcp = Vector3::Cross(ca, cp);
    tempB = Vector3::Dot(triNormal, caXcp) / whole;
    if (tempB < 0.0) {
        return false;
    }

    tempG = 1 - tempA - tempB;
    if (tempG < 0.0) {
        return false;
    }

    // reach here: point is in triangle
    alpha = tempA;
    beta = tempB;
    gamma = tempG;
    return true;
}

bool IntersectedTri(Ray ray, double& t, int& triIndex, double& alpha, double& beta, double& gamma)
{
    bool isIntersected = false;
    for (int i = 0; i < num_triangles; i++) {
        // three vertices of the triangle
        Vector3 a(
            triangles[i].v[0].position[0],
            triangles[i].v[0].position[1],
            triangles[i].v[0].position[2]
        );
        Vector3 b(
            triangles[i].v[1].position[0],
            triangles[i].v[1].position[1],
            triangles[i].v[1].position[2]
        );
        Vector3 c(
            triangles[i].v[2].position[0],
            triangles[i].v[2].position[1],
            triangles[i].v[2].position[2]
        );
        Vector3 triNormal = GetTriangleNormal(a, b, c);
        double nd = Vector3::Dot(triNormal, ray.dir);
        if (Math::NearZero(nd)) {  // no intersection, ray parrallel to plane
            continue;
        }
        else {
            double tMin = Vector3::Dot(triNormal,a) / nd;
            if (tMin <= 0) {
                continue;
            }
            else {
                // closer intersection
                if (tMin < t) {
                    // point-in-triangle testing
                    if (isPointInTri(ray, tMin, alpha, beta, gamma, a, b, c)) {
                        isIntersected = true;
                        t = tMin;
                        triIndex = i;
                    }
                }
            }
        }
    }
    return isIntersected;
}

bool IsInShadow(Ray shadowRay)
{
    shadowRay.dir.Normalize();

    // check sphere intersection
    for (int i = 0; i < num_spheres; i++) {
        Vector3 sphere(
            spheres[i].position[0],
            spheres[i].position[1],
            spheres[i].position[2]
        );
        Vector3 x = shadowRay.startPos - sphere;
        double b = 2.0 * Vector3::Dot(shadowRay.dir, x);
        double c = Vector3::Dot(x, x) - pow(spheres[i].radius,2.0);
        double bSqaureMinus4c = pow(b,2.0) - 4 * c;
        if (bSqaureMinus4c < 0.0) { // abort if negative
            continue;
        }

        double t = (-b - sqrt(bSqaureMinus4c)) / 2.0;
        if (t > 0.0) {
            return true;
        }
    }

    // check tri intersection
    for (int i = 0; i < num_triangles; i++) {
        // three vertices of the triangle
        Vector3 a(
            triangles[i].v[0].position[0],
            triangles[i].v[0].position[1],
            triangles[i].v[0].position[2]
        );
        Vector3 b(
            triangles[i].v[1].position[0],
            triangles[i].v[1].position[1],
            triangles[i].v[1].position[2]
        );
        Vector3 c(
            triangles[i].v[2].position[0],
            triangles[i].v[2].position[1],
            triangles[i].v[2].position[2]
        );
        Vector3 triNormal = GetTriangleNormal(a, b, c);
        triNormal.Normalize();
        double nd = Vector3::Dot(triNormal, shadowRay.dir);
        if (nd <= 0.0) {  // no intersection
            continue;
        }
        else {
            double tMin = Vector3::Dot(triNormal, a - shadowRay.startPos) / nd;
            if (tMin < Math::EPI) {
                continue;
            }
            else {
                double tempA, tempB, tempG;
                // point-in-triangle testing
                if (isPointInTri(shadowRay, tMin, tempA, tempB, tempG, a, b, c)) {
                    return true;
                }
            }
        }
    }

    // reach here: not in the shadow
    return false;
}

Vector3 GetColorSphere(Ray ray, Ray shadowRay, int lightIndex, int index)
{
    Vector3 finalColor(0.0, 0.0, 0.0);

    Vector3 sphereDiffuse(
        spheres[index].color_diffuse[0],
        spheres[index].color_diffuse[1],
        spheres[index].color_diffuse[2]
    );

    Vector3 sphereCenter(
        spheres[index].position[0],
        spheres[index].position[1],
        spheres[index].position[2]
    );
    Vector3 l = shadowRay.dir;
    Vector3 n = shadowRay.startPos - sphereCenter;
    n.Normalize();
    double LdotN = Vector3::Dot(l,n);
    if (LdotN < 0.0) {
        LdotN = 0.0;
    }
    sphereDiffuse *= LdotN;
    finalColor += sphereDiffuse;

    Vector3 sphereSpecular(
        spheres[index].color_specular[0],
        spheres[index].color_specular[1],
        spheres[index].color_specular[2]
    );
    Vector3 r = Vector3::Reflect(l, n);
    r.Normalize();
    Vector3 v = ray.dir * -1.0;
    double RdotV = Vector3::Dot(r, v);
    if (RdotV < 0.0) {
        RdotV = 0.0;
    }
    double sphereShine = spheres[index].shininess;
    double RdotVpowShine = pow(RdotV, sphereShine);
    sphereSpecular *= RdotVpowShine;
    finalColor += sphereSpecular;

    Vector3 lightCol(
        lights[lightIndex].color[0],
        lights[lightIndex].color[1],
        lights[lightIndex].color[2]
    );
    finalColor.x *= lightCol.x;
    finalColor.y *= lightCol.y;
    finalColor.z *= lightCol.z;

    return finalColor;
}

Vector3 GetColorTri(Ray ray, Ray shadowRay, int lightIndex, double alpha, double beta, double gamma, int index)
{
    Vector3 finalColor(0.0, 0.0, 0.0);

    Vector3 triDiffuse(
        triangles[index].v[0].color_diffuse[0] * alpha + triangles[index].v[1].color_diffuse[0] * beta + triangles[index].v[2].color_diffuse[0] * gamma,
        triangles[index].v[0].color_diffuse[1] * alpha + triangles[index].v[1].color_diffuse[1] * beta + triangles[index].v[2].color_diffuse[1] * gamma,
        triangles[index].v[0].color_diffuse[2] * alpha + triangles[index].v[1].color_diffuse[2] * beta + triangles[index].v[2].color_diffuse[2] * gamma
    );
    Vector3 l = shadowRay.dir;
    Vector3 n(
        triangles[index].v[0].normal[0] * alpha + triangles[index].v[1].normal[0] * beta + triangles[index].v[2].normal[0] * gamma,
        triangles[index].v[0].normal[1] * alpha + triangles[index].v[1].normal[1] * beta + triangles[index].v[2].normal[1] * gamma,
        triangles[index].v[0].normal[2] * alpha + triangles[index].v[1].normal[2] * beta + triangles[index].v[2].normal[2] * gamma
    );
    n.Normalize();
    double LdotN = Vector3::Dot(l, n);
    if (LdotN < 0.0) {
        LdotN = 0.0;
    }
    triDiffuse *= LdotN;
    finalColor += triDiffuse;

    Vector3 triSpecular(
        triangles[index].v[0].color_specular[0] * alpha + triangles[index].v[1].color_specular[0] * beta + triangles[index].v[2].color_specular[0] * gamma,
        triangles[index].v[0].color_specular[1] * alpha + triangles[index].v[1].color_specular[1] * beta + triangles[index].v[2].color_specular[1] * gamma,
        triangles[index].v[0].color_specular[2] * alpha + triangles[index].v[1].color_specular[2] * beta + triangles[index].v[2].color_specular[2] * gamma
    );
    Vector3 r = Vector3::Reflect(l, n);
    r.Normalize();
    Vector3 v = ray.dir * -1.0;
    double RdotV = Vector3::Dot(r, v);
    if (RdotV < 0.0) {
        RdotV = 0.0;
    }
    double triShineVal = triangles[index].v[0].shininess * alpha + triangles[index].v[1].shininess * beta + triangles[index].v[2].shininess * gamma;
    Vector3 triShine(triShineVal, triShineVal, triShineVal);
    double triShineFinal = triShine.Length();
    double RdotVpowShine = pow(RdotV, triShineFinal);
    triSpecular *= RdotVpowShine;
    finalColor += triSpecular;

    Vector3 lightCol(
        lights[lightIndex].color[0],
        lights[lightIndex].color[1],
        lights[lightIndex].color[2]
    );
    finalColor.x *= lightCol.x;
    finalColor.y *= lightCol.y;
    finalColor.z *= lightCol.z;

    return finalColor;
}

Vector3 GetPixelColor(Ray ray)
{
    Vector3 color;
    // check sphere intersection
    double tMinSphere = Math::MAX;
    int sphereIndex = -1;
    bool intersectSphere = IntersectedSphere(ray, tMinSphere, sphereIndex);

    // check triangle intersection
    double tMinTri = Math::MAX;
    int triIndex = -1;
    double alpha, beta, gamma;
    bool intersectTri = IntersectedTri(ray, tMinTri, triIndex, alpha, beta, gamma);
    // ray does not hit anything
    if (intersectSphere == false && intersectTri == false) {
        color = Vector3(1.0, 1.0, 1.0); // white background
    }
    // hit sphere or tri or both
    else {
        // determine the closest intersection
        int cloestIntersectionType = -1; // 0 is sphere, 1 is triangle
        double tMin = -1;
        int index = -1;
        if (intersectSphere == true && intersectTri == false) {
            cloestIntersectionType = 0;
            tMin = tMinSphere;
            index = sphereIndex;
        }
        else if (intersectSphere == false && intersectTri == true) {
            cloestIntersectionType = 1;
            tMin = tMinTri;
            index = triIndex;
        }
        else if (intersectSphere == true && intersectTri == true) {
            if (tMinSphere < tMinTri) {
                cloestIntersectionType = 0;
                tMin = tMinSphere;
                index = sphereIndex;
            }
            else if (tMinTri < tMinSphere) {
                cloestIntersectionType = 1;
                tMin = tMinTri;
                index = triIndex;
            }
        }

        // for each light
        for (int i = 0; i < num_lights; i++) {
            // get shadowRay
            Ray shadowRay;
            Vector3 light(
                lights[i].position[0],
                lights[i].position[1],
                lights[i].position[2]
            );
            shadowRay.startPos = ray.dir * tMin;
            shadowRay.dir = light - shadowRay.startPos;
            shadowRay.dir.Normalize();
            // check if in the shadow
            if (IsInShadow(shadowRay)) {
                color += Vector3(0.0, 0.0, 0.0); // black if in shadow
            }
            else {
                if (cloestIntersectionType == 0) { // Sphere Phong shading
                    color += GetColorSphere(ray, shadowRay, i, index);
                }
                else if (cloestIntersectionType == 1) { // Triangle Phong shading
                    color += GetColorTri(ray, shadowRay, i, alpha, beta, gamma, index);
                }
            }
        }
    }

    // add global ambient light
    Vector3 globalAmbient(
        ambient_light[0], ambient_light[1], ambient_light[2]
    );
    color += globalAmbient;
    // clamping
    color.x = Math::Clamp(color.x, 0.0, 1.0);
    color.y = Math::Clamp(color.y, 0.0, 1.0);
    color.z = Math::Clamp(color.z, 0.0, 1.0);
    // to RGB
    color *= 255.0;
    return color;
}

Ray GetRay(double x, double y)
{
    Ray ray;
    // shoot ray from camera at 0,0,0
    ray.startPos.Set(0.0, 0.0, 0.0);
    // calculate the shoot direction of the ray
    double angle = fov / 2.0 * Math::Pi / 180.0;
    double ratio = (double)WIDTH / (double)HEIGHT;
    double planeOffset = 2.0 / HEIGHT;
    double xOffSet = planeOffset * x - ratio;
    double yOffSet = planeOffset * y - 1.0;
    ray.dir.Set(
        tan(angle) * xOffSet,
        tan(angle) * yOffSet,
        -1.0 // -z
    );
    ray.dir.Normalize();
    return ray;
}

void draw_scene()
{
    // for every pixel
    for(unsigned int x=0; x<WIDTH; x++)
    {
        glPointSize(2.0);  
        glBegin(GL_POINTS);
        for(unsigned int y=0; y<HEIGHT; y++)
        {
            // get the ray's direction from camera to that pixel
            Ray ray = GetRay((double)x, double(y));
            // get the color of the pixel
            Vector3 color = GetPixelColor(ray);
            plot_pixel(x, y, color.x, color.y, color.z);
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

