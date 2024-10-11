#include <iostream>
#include <cmath>
#include <fstream>
#include <optional>
#include <utility>
#include <vector>

using namespace std;

class Vector3D {
    public:

    Vector3D(double x, double y, double z) {
        v[0] = x;
        v[1] = y;
        v[2] = z;
    }

    Vector3D() : v{0,0,0} { }

    double getX(){
        return v[0];
    }

    double getY(){
        return v[1];
    }

    double getZ(){
        return v[2];
    }

    friend ostream& operator <<(ostream& out, const Vector3D& obj){
        out << "<" << obj.v[0] << "," << obj.v[1] << "," << obj.v[2] << ">";

        return out;
    }


    friend Vector3D operator *(double scalar, Vector3D obj){
        Vector3D tmp;

        tmp.v[0] = scalar*obj.v[0];
        tmp.v[1] = scalar*obj.v[1];
        tmp.v[2] = scalar*obj.v[2];

        return tmp;
    }

    private:

    double v[3]{};
};

class Sphere {
public:
    Sphere(Vector3D center, double radius, Vector3D color) {
        this->center = center;
        this->radius = radius;
        this->color = color;
    }

    Vector3D getCenter() {
        return center;
    }

    double getRadius() {
        return radius;
    }

    Vector3D getColor() {
        return color;
    }

private:
    Vector3D center;
    double radius;
    Vector3D color;
};

class Light {

public:
    Light(double intensity, string type, Vector3D vector) {
        this->intensity = intensity;
        this->type = type;
        this->vec = vector;
    }

    double getIntensity() {return intensity;}
    Vector3D getVector() {return vec;}
    string getType() {return type;}

private:
    double intensity;
    Vector3D vec;
    string type;
};

// class light_point {
// public:
//     light_point(double intensity, Vector3D position) {
//         this->intensity = intensity;
//         this->position = position;
//         this->type = "point";
//
//     }
//     double getIntensity() {return intensity;}
//     string getType() {return type;}
//     Vector3D getPosition() {return position;}
//
//     private:
//     double intensity;
//     Vector3D position;
//     string type;
// };
//
// class light_directional {
//   public:
//     light_directional(double intensity, Vector3D direction) {
//         this->intensity = intensity;
//         this->direction = direction;
//         this->type = "directional";
//     }
//
//     double getIntensity() {return intensity;}
//     string getType() {return type;}
//     Vector3D getDirection() {return direction;}
//
// private:
//     double intensity;
//     Vector3D direction;
//     string type;
// };

int Cw = 1000;
int Ch = 1000;
int Vw = 1;
int Vh = 1;
int inf = 999999;
double d = 1;
Vector3D origin(0.0, 0.0, 0.0);
Vector3D bColor(0.0, 0.0, 0.0);

vector<Sphere> spheres;
vector<Light> lights;

void imageGen();
Vector3D canvasToViewport(double x,double y);
Vector3D traceRay(Vector3D,Vector3D,double,double);
pair<double, double> InterceptRaySphere(Vector3D, Vector3D, Sphere);
double computeLighting(Vector3D P, Vector3D N);
bool inRange(double min, double max, double value);
double length(Vector3D vec);
void addValues();
int max(double val);

Vector3D vectorAdd(Vector3D a, Vector3D b);

Vector3D scalarMultiply(Vector3D a, double b);
Vector3D scalarDivide(Vector3D a, double b);


int main(int argc, char const *argv[]) {
    addValues();
    imageGen();

    return 0;
}

Vector3D tupleSubtract(Vector3D a, Vector3D b) {
    return {a.getX()-b.getX(),a.getY()-b.getY(),a.getZ()-b.getZ()};
}
double dotProduct(Vector3D a, Vector3D b) {
    return (a.getX()*b.getX())+(a.getY()*b.getY())+(a.getZ()*b.getZ());
}

void imageGen() {

    ofstream ppmFile("raytracer.ppm");

    if (!ppmFile) {
        cerr << "Error opening file " << __FILE__ << endl;
        return;
    }

    ppmFile << "P3\n" << Cw << " " << Ch << "\n255\n";

    for (int y = Cw/2; y > -Cw/2; y--) {
        for (int x = -Ch/2; x < Ch/2; x++) {

            Vector3D D = canvasToViewport(x,y);
            // cout << D.getX() << " " << D.getY() << " " << D.getZ() << ": ";
            Vector3D color = traceRay(origin, D, d, inf);

            int r = max(color.getX());
            int g = max(color.getY());
            int b = max(color.getZ());

            ppmFile << r << " " << g << " " << b << " ";
            //cout << x << " " << y << ": " << r << " " << g << " " << b << " " << endl;
        }
        ppmFile << "\n";  // Newline after each row of pixels
    }

    ppmFile.close();

    std::cout << "PPM file created: output.ppm" << std::endl;
}

void addValues() {
    Sphere red = Sphere(Vector3D(0,-1,3), 1.0, Vector3D(255,0,0));
    Sphere blue = Sphere(Vector3D(2,0,4), 1.0, Vector3D(0,255,0));
    Sphere green = Sphere(Vector3D(-2,0,4), 1.0, Vector3D(0,0,255));
    Sphere yellow = Sphere(Vector3D(0,-5001,0), 5000, Vector3D(255,255,0));
    spheres.push_back(red);
    spheres.push_back(blue);
    spheres.push_back(green);
    spheres.push_back(yellow);

    Light ambient = Light(0.2, "ambient",Vector3D(0,0,0));
    Light point = Light(0.6, "point", Vector3D(2,1,0));
    Light directional = Light(0.2, "directional", Vector3D(1,4,4));
    lights.push_back(ambient);
    lights.push_back(point);
    lights.push_back(directional);

}

Vector3D canvasToViewport(double x, double y) {

    return {x*Vw/Cw, y*Vh/Ch, d};

}

Vector3D traceRay(Vector3D O, Vector3D D, double t_min, double t_max) {
    double closest_t = inf;
    int closest_sphere_index = -1;
    for (int i = 0; i < spheres.size(); i++) {
        pair<double, double> values = InterceptRaySphere(O,D,spheres[i]);
        double t1 = values.first;
        double t2 = values.second;

        if (inRange(t_min, t_max, t1) && t1 < closest_t) {
            closest_t = t1;
            closest_sphere_index = i;
        }

        if (inRange(t_min, t_max, t2) && t2 < closest_t) {
            closest_t = t2;
            closest_sphere_index = i;
        }
    }
    if (closest_sphere_index == -1) {
        return bColor;
    }
    Vector3D P = vectorAdd(scalarMultiply(D, closest_t), O);
    Vector3D N = tupleSubtract(P, spheres[closest_sphere_index].getCenter());
    N = scalarMultiply(N, length(N));
    return scalarMultiply(spheres[closest_sphere_index].getColor(), computeLighting(P, N));

}

pair<double, double> InterceptRaySphere(Vector3D origin, Vector3D D, Sphere sphere) {

    double r = sphere.getRadius();
    Vector3D CO = tupleSubtract(D,sphere.getCenter());

    double a = dotProduct(D,D);
    double b = 2 * dotProduct(D,CO);
    double c = dotProduct(CO,CO) - (double)(r*r);

    double dis = b*b - 4*a*c;

    if (dis < 0) {
        return pair<double, double>{inf,inf};
    }

    double t1 = (-b + sqrt(dis)) / (2*a);
    double t2 = (-b - sqrt(dis)) / (2*a);

    return { t1,t2 };
}

double computeLighting(Vector3D P, Vector3D N) {
    double i = 0.0;
    Vector3D L;
    for (auto & light : lights) {
        if (light.getType() == "ambient") {
            i += light.getIntensity();
        }
        else {
            if (light.getType() == "point") {
                L = tupleSubtract(light.getVector(), P);
            } else {
                L = light.getVector();
            }
        }


        double n_dot_l = dotProduct(N,L);
        if (n_dot_l > 0.0) {
            i += light.getIntensity() * n_dot_l/(length(L) * length(N));
        }
    }
    return i;
}

double length(Vector3D vec) {
    return sqrt(vec.getX()*vec.getX()+vec.getY()*vec.getY()+vec.getZ()*vec.getZ());
}


bool inRange(double min, double max, double value) {
    if (value > min && value < max) {
        return true;
    }
    return false;
}

int max(double val) {
    if (val > 256) {
        return 255;
    }
    if (val < 0) {
        return 0;
    }
    return (int)val;
}


Vector3D scalarMultiply(Vector3D a, double b) {
    return {a.getX() * b, a.getY() * b, a.getZ() * b};
}
Vector3D scalarDivide(Vector3D a, double b) {
    return {a.getX() / b, a.getY() / b, a.getZ() / b};
}


Vector3D vectorAdd(Vector3D a, Vector3D b) {
    return {a.getX()+b.getX(), a.getY()+b.getY(), a.getZ()+b.getZ()};
}


//
// Created by Adhel on 10/8/2024.
//
