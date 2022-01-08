/**
 * ray_tracer.cpp
 * CS230 Assignment 2, Winter 2012
 * -------------------------------
 * Implement ray tracer here.
 */

#define SET_RED(P, C)   (P = (((P) & 0x00ffffff) | ((C) << 24)))
#define SET_GREEN(P, C)  (P = (((P) & 0xff00ffff) | ((C) << 16)))
#define SET_BLUE(P, C) (P = (((P) & 0xffff00ff) | ((C) << 8)))

#include "RayTracer.h"

using namespace std;

const double Object::small_t = 1e-6;
//--------------------------------------------------------------------------------
// utility functions
//--------------------------------------------------------------------------------
double sqr(const double x)
{
    return x * x;
}

Pixel Pixel_Color(const Vector_3D<double>& color)
{
    Pixel pixel = 0;
    SET_RED(pixel, (unsigned char)(min(color.x, 1.0) * 255));
    SET_GREEN(pixel, (unsigned char)(min(color.y, 1.0) * 255));
    SET_BLUE(pixel, (unsigned char)(min(color.z, 1.0) * 255));
    return pixel;
}
//--------------------------------------------------------------------------------
// Shader
//--------------------------------------------------------------------------------
Vector_3D<double> Phong_Shader::
Shade_Surface(const Ray& ray, const Object& intersection_object, const Vector_3D<double>& intersection_point, const Vector_3D<double>& same_side_normal) const
{
    //Variable declarations
    Vector_3D<double> color(0.0, 0.0, 0.0);
    Vector_3D<double> d = color_diffuse;
    Vector_3D<double> s = color_specular;
    Vector_3D<double> N = same_side_normal.Normalized();
    vector<Light*> lights = intersection_object.material_shader->world.lights;
    bool shadow = intersection_object.material_shader->world.enable_shadows;
    double dif, sp;


    //If shadows arent enabled
    if (!shadow) {
        //for each light in the world
        for (unsigned int i = 0; i < lights.size(); i++) {
            //Create the light ray and the reflective ray
            Vector_3D<double> l_ray = (lights[i]->position - intersection_point).Normalized();
            Vector_3D<double> r = (N * 2 * N.Dot_Product(N, ray.endpoint) - ray.endpoint).Normalized();

            //Determine diffuse and specualr components
            sp = pow(max(0.0, l_ray.Dot_Product(r, l_ray)), specular_power);
            dif = max(0.0, N.Dot_Product(N, l_ray));

            //Add to color
            color += lights[i]->Emitted_Light(ray) * d * dif + lights[i]->Emitted_Light(ray) * sp * s;
        }
    }
    else {
        //Create a shadow ray
        Ray shadow;
        shadow.endpoint = intersection_point + same_side_normal * 0.1;

        //for each light in the world
        for (unsigned int i = 0; i < lights.size(); i++) {
            //Initialize the shadow direction and determine if there is an intersection
            shadow.direction = (lights[i]->position - intersection_point).Normalized();
            const Object* shadowed = world.Closest_Intersection(shadow);
            //Create the light ray and the reflective ray
            Vector_3D<double> l_ray = (lights[i]->position - intersection_point).Normalized();
            Vector_3D<double> r = (N * 2 * N.Dot_Product(N, ray.endpoint) - ray.endpoint).Normalized();
            if (!shadowed)
            {
                //Determine diffuse and specualr components
                sp = pow(max(0.0, l_ray.Dot_Product(r, l_ray)), specular_power);
                dif = max(0.0, N.Dot_Product(N, l_ray));

                //Add to color
                color += lights[i]->Emitted_Light(ray) * d * dif + lights[i]->Emitted_Light(ray) * sp * s;
            }
            else {
            }
        }

    }
    return color;
}


Vector_3D<double> Reflective_Shader::
Shade_Surface(const Ray& ray, const Object& intersection_object, const Vector_3D<double>& intersection_point, const Vector_3D<double>& same_side_normal) const
{
    // Variables
    Vector_3D<double> color(0.0, 0.0, 0.0);
    Vector_3D<double> d = color_diffuse;
    Vector_3D<double> s = color_specular;
    Vector_3D<double> N = same_side_normal.Normalized();
    vector<Light*> lights = intersection_object.material_shader->world.lights;
    bool shadow = intersection_object.material_shader->world.enable_shadows;
    double dif, sp;
    Ray reflection;

    //set up the reflection ray
    reflection.endpoint = intersection_point + same_side_normal * 0.1;
    reflection.direction = ray.direction - (N * ((N.Dot_Product(N, ray.direction)) * 2));

    //If there arent shadows
    if (!shadow) {
        //for each light source in the world
        for (unsigned int i = 0; i < lights.size(); i++) {
            //creates the light ray and the reflection ray
            Vector_3D<double> l_ray = (lights[i]->position - intersection_point).Normalized();
            Vector_3D<double> r = (N * 2 * N.Dot_Product(N, ray.endpoint) - ray.endpoint).Normalized();

            //calculates the diffuse and specular components
            sp = pow(max(0.0, l_ray.Dot_Product(r, l_ray)), specular_power);
            dif = max(0.0, N.Dot_Product(N, l_ray));

            //if the recursion depth is not higher than the max depth
            if (ray.recursion_depth < intersection_object.material_shader->world.recursion_depth_limit) {
                //Cast another ray from the reflection (RECURSIVE PART)
                color += intersection_object.material_shader->world.Cast_Ray(reflection, ray) * reflectivity;
            }
            color += lights[i]->Emitted_Light(ray) * d * dif + lights[i]->Emitted_Light(ray) * sp * s;

        }
    }
    //if shadows are enabled
    else {

        Ray shadow;
        shadow.endpoint = intersection_point + same_side_normal * 0.1;
        //for each light source in the world
        for (unsigned int i = 0; i < lights.size(); i++) {
            //create the shadow ray, light ray and reflection ray
            shadow.direction = (lights[i]->position - intersection_point).Normalized();
            const Object* shadowed = world.Closest_Intersection(shadow);
            Vector_3D<double> l_ray = (lights[i]->position - intersection_point).Normalized();
            Vector_3D<double> r = (N * 2 * N.Dot_Product(N, ray.endpoint) - ray.endpoint).Normalized();

            //if there is an intersection with the shadow
            if (!shadowed)
            {
                //compute specular and diffuse component
                sp = pow(max(0.0, l_ray.Dot_Product(r, l_ray)), specular_power);
                dif = max(0.0, N.Dot_Product(N, l_ray));

                //if the recursion depth is not higher than the max depth
                if (ray.recursion_depth < intersection_object.material_shader->world.recursion_depth_limit) {

                    //Cast another ray from the reflection (RECURSIVE PART)	//Cast another ray from the reflection (RECURSIVE PART)
                    color += intersection_object.material_shader->world.Cast_Ray(reflection, ray) * reflectivity;
                }
                color += lights[i]->Emitted_Light(ray) * d * dif + lights[i]->Emitted_Light(ray) * sp * s;
            }
        }

    }
    return color;
}

Vector_3D<double> Flat_Shader::
Shade_Surface(const Ray& ray, const Object& intersection_object, const Vector_3D<double>& intersection_point, const Vector_3D<double>& same_side_normal) const
{
    return color;
}

//--------------------------------------------------------------------------------
// Objects
//--------------------------------------------------------------------------------
// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Sphere::Intersection(Ray& ray) const
{

    Vector_3D<double> p0 = ray.direction;
    Vector_3D<double> p1 = ray.endpoint;

    double a = p0.Dot_Product(p0, p0);
    double b = 2 * p0.Dot_Product(p0, p1 - center);
    double c = center.Dot_Product(p1 - center, p1 - center) - sqr(radius);

    double desc = b * b - 4 * a * c;

    //No Intersection
    if (desc < 0) {
        return false;
    }
    // Intersection
    else {

        if ((-b - sqrt(desc)) / (2 * a) < 0 && (-b + sqrt(desc)) / (2 * a) < 0) {
            return false;
        }

        ray.t_max = min((-b - sqrt(desc)) / (2 * a), (-b + sqrt(desc)) / (2 * a));


        ray.current_object = this;
        ray.semi_infinite = false;
        return true;
    }

}

Vector_3D<double> Sphere::Normal(const Vector_3D<double>& location) const
{
    Vector_3D<double> normal;
    normal = location - center;
    return normal.Normalized();
}

// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Plane::Intersection(Ray& ray) const
{
    double t1 = normal.Dot_Product(normal, ray.endpoint - x1) * -1;
    double t2 = normal.Dot_Product(normal, ray.direction);
    //No intersection
    if (t1 == 0 && t2 == 0) {
        return false;
    }
    else if (t1 == 0) {
        return false;
    }
    else if (t2 == 0) {
        return false;
    }
    else {
        if (t1 / t2 < 0) {
            return false;
        }
        //Intersection
        ray.t_max = t1 / t2;
        ray.current_object = this;
        ray.semi_infinite = false;
        return true;
    }

}

Vector_3D<double> Plane::Normal(const Vector_3D<double>& location) const
{
    return normal;
}
//--------------------------------------------------------------------------------
// Camera
//--------------------------------------------------------------------------------
// Find the world position of the input pixel
Vector_3D<double> Camera::World_Position(const Vector_2D<int>& pixel_index)
{
    Vector_3D<double> result;
    Vector_2D<double> a = film.pixel_grid.X(pixel_index);
    result = focal_point + horizontal_vector * a.x + vertical_vector * a.y;
    return result - position;
}
//--------------------------------------------------------------------------------
// Render_World
//--------------------------------------------------------------------------------
// Find the closest object of intersection and return a pointer to it
//   if the ray intersects with an object, then ray.t_max, ray.current_object, and ray.semi_infinite will be set appropriately
//   if there is no intersection do not modify the ray and return 0
const Object* Render_World::
Closest_Intersection(Ray& ray)
{
    double temp = 100000000;
    Object* temp_obj = 0;
    //for all the objects in the world, find which is closest
    for (unsigned int i = 0; i < objects.size(); i++) {
        if (objects[i]->Intersection(ray)) {
            if (ray.t_max < temp) {
                temp = ray.t_max;
                temp_obj = objects[i];
            }
        }

    }
    ray.t_max = temp;
    ray.current_object = temp_obj;
    return temp_obj;
}

// set up the initial view ray and call 
void Render_World::Render_Pixel(const Vector_2D<int>& pixel_index)
{
    //set up the initial view ray
    Ray ray;
    ray.endpoint = camera.position;
    ray.direction = camera.World_Position(pixel_index);
    Ray dummy_root;
    Vector_3D<double> color = Cast_Ray(ray, dummy_root);
    camera.film.Set_Pixel(pixel_index, Pixel_Color(color));
}

// cast ray and return the color of the closest intersected surface point, 
// or the background color if there is no object intersection
Vector_3D<double> Render_World::Cast_Ray(Ray& ray, const Ray& parent_ray)
{
    Vector_3D<double> color;
    const Object* closest = Closest_Intersection(ray);

    ray.recursion_depth++;
    if (closest == NULL) {
        return color;
    }

    color = closest->material_shader->Shade_Surface(ray, *closest, ray.Point(ray.t_max), closest->Normal(ray.Point(ray.t_max)));
    return color;
}
