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


#include "vec3.h"
#include <fstream>
#include <string>
#include <iostream>



struct ray {
    vec3 origin;
    vec3 direction; 
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
            ray r;
            r.origin = origin;
            r.direction = lower_left_corner + s*horizontal + t*vertical - origin;
            return r;
        }

        vec3 origin;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
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
};

struct sphere {
    vec3 center;
    float radius;
    vec3 pigment;
    surface surfac;
};

float eps = 0.001;
float hit(ray r, sphere s) {
    float t_min = 0.001;
    float t_max = MAXFLOAT;
    vec3 center = s.center;
    float radius = s.radius;

    vec3 oc = r.origin - center;
    float a = dot(r.direction, r.direction);
    float b = dot(oc, r.direction);
    float c = dot(oc, oc) - radius*radius;
    float discriminant = b*b - 4*a*c;
    if (discriminant > 0-eps) {
        float temp = (-b - sqrt(discriminant))/a;
        float temp2 = (-b + sqrt(discriminant)) / a;
        if (temp + eps > 0 && temp2 + eps > temp){
            return temp;
        }
        return -1.0;
    }
    return -1.0;;

}
/*
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
*/


int MAX = 4;



light light_sources[20];
vec3 pigments[25];
surface surfaces[25];
sphere objects[100];

int number_of_light = 0;
int number_of_pigments = 0;
int number_of_surface = 0;

int number_of_objects = 0;
vec3 background_color = vec3(0.5,0.5,0.5);



vec3 phong(light l, vec3 point, vec3 normal,sphere s,ray r){

    
    float a = l.attenuation[0];
    float b = l.attenuation[1];
    float c = l.attenuation[2];
    vec3 position = l.position;
   
    float surface_ka = s.surfac.ka;
    float surface_kd = s.surfac.kd;
    float surface_ks = s.surfac.ks;
    float surface_alfa = s.surfac.alfa;
    float surface_kr = s.surfac.kr;
    normal = unit_vector(normal);
    vec3 light_vector = unit_vector(position - point);
    float d = sqrt(dot(light_vector,light_vector));
    float att = 1/(a+b*d+c*d*d) ;
    vec3 ambient_light_color = light_sources[0].rgb;
    vec3 color_a = surface_ka*vec3(s.pigment[0]*ambient_light_color[0],s.pigment[1]*ambient_light_color[1],s.pigment[2]*ambient_light_color[2]);
    float l_n = dot(light_vector,normal);
    if(l_n < 0)
        l_n = 0;
    
    vec3 color_d = (l_n*surface_kd*vec3(s.pigment[0]*l.rgb[0],s.pigment[1]*l.rgb[1],s.pigment[2]*l.rgb[2])) / att;
    vec3 v = unit_vector(r.origin - point);
    vec3 h = unit_vector(light_vector + v);
    float h_n = dot(h,normal);
    if(h_n < 0)
        h_n = 0;
    float i_s = pow(h_n,surface_alfa);
    vec3 color_s = (i_s*surface_ks*vec3(s.pigment[0]*l.rgb[0],s.pigment[1]*l.rgb[1],s.pigment[2]*l.rgb[2])) / att ;
    return color_a + color_d + color_s;
}


bool visible(vec3 P,light l){
    //dot product of normal and light
    vec3 position = l.position;
    if (dot(P,position) + eps > 0);  
        return 1;
    return 0;
    
}
vec3 color(ray r, int depth) {
    float t = MAXFLOAT;
    int closest_object_index = -1;
    for(int i = 0; i < number_of_objects; i ++){
        sphere s = objects[i];
        float a = hit(r,s);
        if(t>a){
            t = a;
            closest_object_index = i;
        }
    }

    vec3 localC = vec3(0.0,0.0,0.0), reflectedC, transmittedC;
    vec3 P; vec3 normal;
    P = r.origin + t*r.direction;
    sphere s = objects[closest_object_index];
    normal = r.origin + t*r.direction-s.center;
    
    for(int i = 1; i < number_of_light; i++){
        light l = light_sources[i];
        if (visible(P,l)){
            localC += phong(l, P, normal,s,r);
            return localC;
        }
        // add reflectivituy
    }
}

/*

light light_sources[20];
vec3 pigments[25];
surface surfaces[25];
sphere objects[100];

int number_of_light = 0;
int number_of_pigments = 0;
int number_of_surface = 0;

int number_of_objects = 0;
vec3 background_color = vec3(0.5,0.5,0.5);

 */

int nx = 1000;
int ny = 1000;
vec3 lookfrom(20,9,0);
vec3 lookat(0,0,-3);
int aspect = 90;

void init_scene() {
    std::cout << "Hello World!\n";
    std::cout << "I am learning C++\n";
    
    number_of_light = 3;
    /*for(int i=0;i < number_of_light; i++){
        light l;
        l.position = vec3(0,0,0);
        l.rgb = vec3(0,0,0);
        l.attenuation = vec3(0,0,0);
        light_sources[i] = l;
    }
     */
    light l;
    l.position = vec3(0,0,0);
    l.rgb = vec3(1,1,1);
    l.attenuation = vec3(1,0,0);
    light_sources[0] = l;

    l.position = vec3(0,0,0);
    l.rgb = vec3(1,1,1);
    l.attenuation = vec3(1,0,0);
    light_sources[1] = l;

    l.position = vec3(0,0,0);
    l.rgb = vec3(1,1,1);
    l.attenuation = vec3(1,0,0);
    light_sources[2] = l;
    
    number_of_pigments = 3;
    /*for(int i=0; i<number_of_pigments; i++){
        vec3 pigm(0.0,0.0,0.0);
        pigments[i] = pigm;
    } */
    pigments[0] = vec3(1.0,0.0,0.0);
    pigments[1] = vec3(0.0,1.0,0.0);
    pigments[2] = vec3(0.0,0.0,1.0);


    /*
    number_of_surface = 2;
    for(int i=0; i<number_of_surface; i++){
        float ka;
        float kd;
        float ks;
        float alfa;
        float kr;
        surface s;
        s.ka = ka;
        s.kd = kd;
        s.ks = ks;
        s.alfa = alfa;
        s.kr = kr;
        
        surfaces[i] = s;
    } */

    float ka = 0.4;
    float kd = 0.6;
    float ks = 0.0;
    float alfa = 1.0;
    float kr = 0.0;
    surface s;
    s.ka = ka;
    s.kd = kd;
    s.ks = ks;
    s.alfa = alfa;
    s.kr = kr;
    surfaces[0] = s;

    ka = 0.4;
    kd = 0.6;
    ks = 0.7;
    alfa = 500.0;
    kr = 0.0;

    s.ka = ka;
    s.kd = kd;
    s.ks = ks;
    s.alfa = alfa;
    s.kr = kr;    
    surfaces[1] = s;

    number_of_objects = 4;
    /*for(int i=0; i<number_of_objects; i++){
        sphere s;
        s.center = vec3();
        s.pigment = pigments[0];
        s.surfac = surfaces[0];
        s.radius = 3.0;
        objects[i];
    } */
    sphere sp;
    sp.center = vec3();
    sp.pigment = pigments[1];
    sp.surfac = surfaces[0];
    sp.radius = 3.0;
    objects[0];
    
    sp.center = vec3();
    sp.pigment = pigments[0];
    sp.surfac = surfaces[1];
    sp.radius = 3.0;
    objects[1];

    sp.center = vec3();
    sp.pigment = pigments[2];
    sp.surfac = surfaces[1];
    sp.radius = 3.0;
    objects[2];

    sp.center = vec3();
    sp.pigment = pigments[2];
    sp.surfac = surfaces[1];
    sp.radius = 3.0;
    objects[3];
    std::cout << "I am learning C++\n";

}


int main() {
    
    int ns = 10;
    //std::cout << "P3\n" << nx << " " << ny << "\n255\n";
    init_scene();

    camera cam(lookfrom, lookat, vec3(0,1,0), aspect, float(nx)/float(ny));

    for (int j = ny-1; j >= 0; j--) {
        for (int i = 0; i < nx; i++) {
            vec3 col(0, 0, 0);
            for (int s=0; s < ns; s++) {
                float u = float(i) / float(nx);
                float v = float(j) / float(ny);
                ray r = cam.get_ray(u, v);
                col += color(r,0);
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
