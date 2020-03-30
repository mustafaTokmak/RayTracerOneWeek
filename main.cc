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
#include <string.h>


int debug = 0;

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

    vec3 oc = center - r.origin ;
    float a = dot(r.direction, r.direction);
    float b = (-2)*dot(oc, r.direction);
    float c = dot(oc, oc) - radius*radius;
    float discriminant = b*b - 4*a*c;
    if(debug){
        std::cout << "discriminant" << discriminant << "\n "  ;
    }
    if (discriminant > 0-eps) {
        float temp = (-b - sqrt(discriminant))/a;
        float temp2 = (-b + sqrt(discriminant)) / a;
        if(debug){
            std::cout << "temp " << temp << "\n "  ;
            std::cout << "temp2 " << temp2 << "\n "  ;
        }
        if (temp + eps > 0 ){
            return temp;
        }
        else if(temp2 + eps > 0){
            return temp2;
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
    if(debug){
        std::cout << "normal  " << normal << "\n "  ;
    }
    vec3 light_vector = (position - point);
    if(debug){
        std::cout << "light_vector  " << light_vector << "\n "  ;
    }

    float d = sqrt(dot(light_vector,light_vector));
    if(debug){
        std::cout << "d  " << d << "\n "  ;
    }

    float att = 1.0/(a+b*d+c*d*d) ;
    if(debug){
        std::cout << "att  " << att << "\n "  ;
    }
    
    vec3 ambient_light_color = light_sources[0].rgb;
    if(debug){
        std::cout << "ambient_light_color  " << ambient_light_color << "\n "  ;
    }
    if(debug){
        std::cout << "s.pigment  " << s.pigment << "\n "  ;
    }
    
    vec3 color_a = surface_ka*vec3(s.pigment[0]*ambient_light_color[0],s.pigment[1]*ambient_light_color[1],s.pigment[2]*ambient_light_color[2]);
    if(debug){
        std::cout << "color_a  " << color_a << "\n "  ; 
    }
    light_vector = unit_vector(light_vector);
    if(debug){
        std::cout << "light_vector  " << light_vector << "\n "  ;
    }
    
    float l_n = dot(light_vector,normal);
    if(debug){
        std::cout << "l_n  " << l_n << "\n "  ;
    }
    if(l_n < 0)
        l_n = 0;
    if(debug){
        std::cout << "l_n  " << l_n << "\n "  ;
    }
    vec3 color_d = (l_n*surface_kd*vec3(s.pigment[0]*l.rgb[0],s.pigment[1]*l.rgb[1],s.pigment[2]*l.rgb[2])) / att;
    if(debug){
        std::cout << "color_d  " << color_d << "\n "  ;
    }
    vec3 v = unit_vector(r.origin - point);
    if(debug){
        std::cout << "v  " << v << "\n "  ;
    }
    vec3 h = unit_vector(light_vector + v);
    if(debug){
        std::cout << "h  " << h << "\n "  ;
    }
    
    float h_n = dot(h,normal);
    if(h_n < 0)
        h_n = 0;
    if(debug){
        std::cout << "h_n  " << h_n << "\n "  ;
    }
    float i_s = pow(h_n,surface_alfa);
    if(debug){
        std::cout << "i_s  " << i_s << "\n "  ;
    }
    vec3 color_s = (i_s*surface_ks*vec3(s.pigment[0]*l.rgb[0],s.pigment[1]*l.rgb[1],s.pigment[2]*l.rgb[2])) / att ;
    if(debug){
        std::cout << "color_s  " << color_s << "\n "  ;
    }
    if(debug){
        std::cout << "color_a + color_d + color_s  " << color_a + color_d + color_s << "\n "  ;
    }
    return color_a + color_d + color_s;
}


int visible(vec3 P,light l,vec3 normal,int index){
    //dot product of normal and light
    
    normal = unit_vector(normal);
    vec3 position = l.position;
    vec3 light_vector = unit_vector(position-P);
    float a  = dot(light_vector,normal)-eps;
    if(debug){
        std::cout << "visible_a " << a <<"  " << index << "\n "  ;
    }
    if ( a  < 0) {
        return -1;
    }
    ray r;

    r.origin = P;
    r.direction = light_vector;
    for(int i = 0; i < number_of_objects; i ++){
        if(i == index){
            continue;
        }
        sphere s = objects[i];
        float a = hit(r,s);
        
        if(a > 0 ){
            return -1;
        }
    } 
    return 1;
    
}
vec3 color(ray r, int depth) {
    float t = MAXFLOAT;
    int closest_object_index = -1;
    for(int i = 0; i < number_of_objects; i ++){
        sphere s = objects[i];
        float a = hit(r,s);
        if(debug && a > 0){
            std::cout << "hit " << a << "\n "  ;
        }
        if(a > 0 && t>a){
            t = a;
            closest_object_index = i;
        }
        if(debug && a > 0 ){
            std::cout << i << "\n"  ;
        }
    }
    if(closest_object_index == -1){
        return background_color;
    }
    
    vec3 localC = vec3(0.0,0.0,0.0), reflectedC, transmittedC;
    vec3 P; vec3 normal;
    P = r.origin + t*r.direction;
    sphere s = objects[closest_object_index];
    normal = r.origin + t*r.direction-s.center;
    if(debug){
        std::cout << "normal " <<normal << "\n "  ;
    }
    //ambient
    
    for(int i = 1; i < number_of_light; i++){
        
        light l = light_sources[i];
        if (visible(P,l,normal,closest_object_index) >0){
            localC += phong(l, P, normal,s,r);
            
        }
        else{
            l = light_sources[0];
            localC += phong(l, P, normal,s,r);
        }
        // add reflectivituy
    }
    return localC;
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
vec3 lookfrom(0.0,0.0,0.0);
vec3 lookat(0,0,-1);
vec3 up(0.0,1.0,0.0);
int fovy = 90;

void init_scene() {
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

    l.position = vec3(10,100,10);
    l.rgb = vec3(1,1,1);
    l.attenuation = vec3(1,0,0);
    light_sources[1] = l;

    l.position = vec3(100,100,100);
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
    if(debug){
                std::cout << "pigments" << "\n" << pigments[0] << "\n "  ;
        }

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
    sp.center = vec3(1,0,-8);
    
    sp.pigment = pigments[1];
    sp.surfac = surfaces[0];
    sp.radius = 2;
    objects[0] = sp;
    
    sp.center = vec3(3,5,-10);
    sp.pigment = pigments[0];
    sp.surfac = surfaces[1];
    sp.radius = 3.0;
    objects[1] = sp;

    sp.center = vec3(10,-5,-25);
    sp.pigment = pigments[2];
    sp.surfac = surfaces[1];
    sp.radius = 10.0;
    objects[2] = sp;

    sp.center = vec3(-10,0,-25);
    sp.pigment = pigments[2];
    sp.surfac = surfaces[1];
    sp.radius = 10.0;
    objects[3]= sp;

    
    

     
    
}

int main() {
    debug = 1;
    int ns = 10;
    //std::cout << "P3\n" << nx << " " << ny << "\n255\n";
    //init();
    init_scene();

    camera cam(lookfrom, lookat, up, fovy, float(nx)/float(ny));

    for (int j = ny-1; j >= 0; j--) {
        for (int i = 0; i < nx; i++) {
            vec3 col(0, 0, 0);
            
            float u = float(i) / float(nx);
            float v = float(j) / float(ny);
            if(debug){
                u = float(500) / float(nx);
                v = float(500) / float(ny);
            }
            
            ray r = cam.get_ray(u, v);
            if(debug){
                std::cout << r.direction << "\n" << r.origin << "\n "  ;
            }
                
            col = color(r,0);
            //col /= float(ns);
            //col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
            
            int ir = int(255.99*col[0]);
            int ig = int(255.99*col[1]);
            int ib = int(255.99*col[2]);
            std::cout << ir << " " << ig << " " << ib << "\n";
            if(debug){
                exit(0);
            }
        }
    }
}
