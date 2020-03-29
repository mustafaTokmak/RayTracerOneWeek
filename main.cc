//==================================================================================================
// Written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is distributed
// without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication along
// with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==================================================================================================

#include <iostream>
#include "hitable_list.h"
#include "float.h"
#include "camera.h"
#include "material.h"



class ray
{
    public:
        ray() {}
        ray(const vec3& a, const vec3& b) { A = a; B = b; }
        vec3 origin() const       { return A; }
        vec3 direction() const    { return B; }
        vec3 point_at_parameter(float t) const { return A + t*B; }

        vec3 A;
        vec3 B;
};


class camera {
    public:
        camera(vec3 lookfrom, vec3 lookat, vec3 vup, float vfov, float aspect) {
            // vfov is top to bottom in degrees
            float theta = vfov*M_PI/180;
            float half_height = tan(theta/2);
            float half_width = aspect * half_height;
            origin = lookfrom;
            w = unit_vector(lookfrom - lookat);
            u = unit_vector(cross(vup, w));
            v = cross(w, u);
            lower_left_corner = origin  - half_width*u -half_height*v - w;
            horizontal = 2*half_width*u;
            vertical = 2*half_height*v;
        }
        ray get_ray(float s, float t) {
            return ray(origin, lower_left_corner + s*horizontal + t*vertical - origin);
        }

        vec3 origin;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
};

struct hit_record
{
    float t;
    vec3 p;
    vec3 normal;
    material *mat_ptr;
};

struct surface
{
    float ka;
    float kd;
    float ks;
    float alfa;
    float kr;
};

struct light
{
    vec3 position;
    vec3 rgb;
    vec3 attenuation;
}



class hitable  {
    public:
        virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const = 0;
};


class sphere: public hitable  {
    public:
        sphere() {}
        sphere(vec3 cen, float r, vec3 pigm,surface s) : center(cen), radius(r), surface(s)  {};
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center;
        float radius;
        float pigment;
        surface surface;
};


bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;
    float a = dot(r.direction(), r.direction());
    float b = dot(oc, r.direction());
    float c = dot(oc, oc) - radius*radius;
    float discriminant = b*b - a*c;
    if (discriminant > 0) {
        float temp = (-b - sqrt(discriminant))/a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            return true;
        }
        temp = (-b + sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            return true;
        }
    }
    return false;
}

vec3 background_color = vec3(0.5,0.5,0.5);

int MAX = 4;



light light_sources[20];
vec3 pigments[25];
surface surfaces[25];


int number_of_light = 0;
int number_of_pigments = 0;
int number_of_surface = 0;

int number_of_objects = 0;



vec3 compute_normal(vec3 P){ 
    return vec3(0.0,0.0,0.0);
}
bool visible(vec3 P,light l){
    return 1;
}

vec3 phong(light l, vec3 point, vec3 normal){
    return background_color;
}

/*vec3 trace(const ray& r,int depth,spehere){
    vec3 localC = vec3(0.0,0.0,0.0), reflectedC, transmittedC;
    vec3 P; vec3 normal;
    normal = compute_normal(P);
    for(int i = 0; i < number_of_light; i++){
        light l = light_sources[i];

        if (visible(P,l)){
            localC += phong(l, P, normal);
            return background_color;
        }
        // add reflectivituy
    }

} */

vec3 color(const ray& r, hitable *world, int depth) {
    


    hit_record rec;
    if (world->hit(r, 0.001, MAXFLOAT, rec)) {
        vec3 localC = vec3(0.0,0.0,0.0), reflectedC, transmittedC;
        vec3 P; vec3 normal;
        normal = compute_normal(P);
        for(int i = 0; i < number_of_light; i++){
            light l = light_sources[i];
            if (visible(P,l)){
                localC += phong(l, P, normal);
                return localC;
            }
            // add reflectivituy
        }


        if (depth < MAX) {
             return background_color;//attenuation*color(scattered, world, depth+1);
        }
        else {
            return background_color;
        }
    }
    else {
        return background_color;
    }
}



class hitable_list: public hitable  {
    public:
        hitable_list() {}
        hitable_list(hitable **l, int n) { list = l; list_size = n; }
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        hitable **list;
        int list_size;
};

bool hitable_list::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
    hit_record temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;
    for (int i = 0; i < list_size; i++) {
        if (list[i]->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    return hit_anything;
}

hitable *random_scene() {
    int n = 4;
    hitable **list = new hitable*[n+1];

    list[0] =  new sphere(vec3(0,-1000,0), 1000, new lambertian(vec3(0.5, 0.5, 0.5)));

    int i = 1;
    list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
    list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
    list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));
    

    return new hitable_list(list,i);
}


int main() {
    int nx = 1200;
    int ny = 800;
    int ns = 10;
    std::cout << "P3\n" << nx << " " << ny << "\n255\n";
    hitable *world = random_scene();

    vec3 lookfrom(20,9,0);
    vec3 lookat(0,0,-3);

    camera cam(lookfrom, lookat, vec3(0,1,0), 60, float(nx)/float(ny));

    for (int j = ny-1; j >= 0; j--) {
        for (int i = 0; i < nx; i++) {
            vec3 col(0, 0, 0);
            for (int s=0; s < ns; s++) {
                float u = float(i + random_double()) / float(nx);
                float v = float(j + random_double()) / float(ny);
                ray r = cam.get_ray(u, v);
                col += color(r, world,0);
            }
            col /= float(ns);
            col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
            int ir = int(255.99*col[0]);
            int ig = int(255.99*col[1]);
            int ib = int(255.99*col[2]);
            std::cout << ir << " " << ig << " " << ib << "\n";
        }
    }
}
