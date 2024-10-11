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
    friend Vector3D operator *(Vector3D obj, double scalar){
        Vector3D tmp;

        tmp.v[0] = obj.v[0]*scalar;
        tmp.v[1] = obj.v[1]*scalar;
        tmp.v[2] = obj.v[2]*scalar;

        return tmp;
    }
    friend Vector3D operator *(Vector3D a, Vector3D b) {
        Vector3D tmp;

        tmp.v[0] = a.v[0]*b.v[0];
        tmp.v[1] = a.v[1]*b.v[1];
        tmp.v[2] = a.v[2]*b.v[2];

        return tmp;
    }

    friend Vector3D operator /(double scalar, Vector3D obj) {
        Vector3D tmp;

        tmp.v[0] = scalar/obj.v[0];
        tmp.v[1] = scalar/obj.v[1];
        tmp.v[2] = scalar/obj.v[2];

        return tmp;
    }
    friend Vector3D operator /(Vector3D obj, double scalar) {
        Vector3D tmp;

        tmp.v[0] = obj.v[0]/scalar;
        tmp.v[1] = obj.v[1]/scalar;
        tmp.v[2] = obj.v[2]/scalar;

        return tmp;
    }
    friend Vector3D operator /(Vector3D a, Vector3D b) {
        Vector3D tmp;

        tmp.v[0] = a.v[0]/b.v[0];
        tmp.v[1] = a.v[1]/b.v[1];
        tmp.v[2] = a.v[2]/b.v[2];

        return tmp;
    }

    friend Vector3D operator -(double scalar, Vector3D obj) {
        Vector3D tmp;

        tmp.v[0] = obj.v[0]-scalar;
        tmp.v[1] = obj.v[1]-scalar;
        tmp.v[2] = obj.v[2]-scalar;

        return tmp;
    }
    friend Vector3D operator -(Vector3D obj, double scalar) {
        Vector3D tmp;

        tmp.v[0] = scalar-obj.v[0];
        tmp.v[1] = scalar-obj.v[1];
        tmp.v[2] = scalar-obj.v[2];

        return tmp;
    }
    friend Vector3D operator -(Vector3D a, Vector3D b) {
        Vector3D tmp;

        tmp.v[0] = a.v[0]-b.v[0];
        tmp.v[1] = a.v[1]-b.v[1];
        tmp.v[2] = a.v[2]-b.v[2];

        return tmp;
    }

    friend Vector3D operator +(double scalar, Vector3D obj) {
        Vector3D tmp;

        tmp.v[0] = obj.v[0]+scalar;
        tmp.v[1] = obj.v[1]+scalar;
        tmp.v[2] = obj.v[2]+scalar;

        return tmp;
    }
    friend Vector3D operator +(Vector3D obj, double scalar) {
        Vector3D tmp;

        tmp.v[0] = scalar+obj.v[0];
        tmp.v[1] = scalar+obj.v[1];
        tmp.v[2] = scalar+obj.v[2];

        return tmp;
    }
    friend Vector3D operator +(Vector3D a, Vector3D b) {
        Vector3D tmp;

        tmp.v[0] = a.v[0]+b.v[0];
        tmp.v[1] = a.v[1]+b.v[1];
        tmp.v[2] = a.v[2]+b.v[2];

        return tmp;
    }

    private:

    double v[3]{};
};

class Sphere {
public:
    Sphere(Vector3D center, double radius, Vector3D color, double spec) {
        this->center = center;
        this->radius = radius;
        this->color = color;
        this->spec = spec;
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
    double getSpec() {
        return spec;
    }

private:
    Vector3D center;
    double radius;
    Vector3D color;
    double spec;
};

class Light {

public:
    Light(double intensity, string type, Vector3D vec) {
        this->intensity = intensity;
        this->type = std::move(type);
        this->vec = vec;
    }

    double getIntensity() {return intensity;}
    Vector3D getVector() {return vec;}
    string getType() {return type;}

private:
    double intensity;
    Vector3D vec;
    string type;
};

int Cw = 300;
int Ch = 300;
int Vw = 1;
int Vh = 1;
int inf = 9999999;
double d = 1;
Vector3D origin(0.0, 0.0, 0.0);
Vector3D bColor(0.0, 0.0, 0.0);

vector<Sphere> spheres;
vector<Light> lights;

void imageGen();
Vector3D canvasToViewport(double x,double y);
Vector3D traceRay(Vector3D,Vector3D,double,double);
pair<double, double> InterceptRaySphere(Vector3D, Vector3D, Sphere);
double computeLighting(Vector3D P, Vector3D N, Vector3D V, double s);
bool inRange(double min, double max, double value);
double length(Vector3D vec);
void addValues();
int max(double val);

Vector3D vectorAdd(Vector3D a, Vector3D b);
Vector3D vectorSub(Vector3D a, Vector3D b);

Vector3D scalarMultiply(Vector3D a, double b);
Vector3D scalarDivide(Vector3D a, double b);


int main(int argc, char const *argv[]) {
    addValues();
    imageGen();


    // Vector3D a(1, 1, 1);
    // Vector3D b(10, 10, 10);
    // double s = 2;
    //
    // cout << "a: " << a << endl;
    // cout << "b: " << b << endl;
    // cout << "s: " << s << endl;
    //
    // Vector3D c = a-b;
    //
    // cout << c << endl;

    return 0;
}


double dotProduct(Vector3D a, Vector3D b) {
    return a.getX()*b.getX()+a.getY()*b.getY()+a.getZ()*b.getZ();
}

void imageGen() {

    ofstream ppmFile("raytracer.ppm");

    if (!ppmFile) {
        cerr << "Error opening file " << __FILE__ << endl;
        return;
    }

    ppmFile << "P3\n" << Cw << " " << Ch << "\n255\n";

    for (int x = -Cw/2; x < Cw/2; x++) {
        for (int y = -Ch/2; y < Ch/2; y++) {

            Vector3D D = canvasToViewport(y,-x);
            // cout << D.getX() << " " << D.getY() << " " << D.getZ() << ": ";
            Vector3D color = traceRay(origin, D, d, inf);

            int r = max(color.getX());
            int g = max(color.getY());
            int b = max(color.getZ());

            ppmFile << r << " " << g << " " << b << " ";
            // cout << x << " " << y << ": " << r << " " << g << " " << b << " " << endl;
        }
        ppmFile << "\n";  // Newline after each row of pixels
    }

    ppmFile.close();

    cout << "PPM file created: output.ppm" << endl;
}

void addValues() {
    Sphere red = Sphere(Vector3D(0,-1,3), 1.0, Vector3D(255,0,0), 500);
    Sphere blue = Sphere(Vector3D(-2,0,4), 1.0, Vector3D(0,255,0), 10);
    Sphere green = Sphere(Vector3D(2,0,4), 1.0, Vector3D(0,0,255), 500);
    Sphere yellow = Sphere(Vector3D(0,-5001,0), 5000, Vector3D(255,255,0), 1000);
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
    // cout << spheres.size();
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
    //Vector3D P = vectorAdd(scalarMultiply(D, closest_t), O);
    Vector3D P = O + closest_t * D;
    Vector3D N = P - spheres[closest_sphere_index].getCenter();
    //Vector3D N = vectorSub(P, spheres[closest_sphere_index].getCenter());
    //N = scalarMultiply(N, length(N));
    N = N/length(N);

    return spheres[closest_sphere_index].getColor() * computeLighting(P, N, -1*D, spheres[closest_sphere_index].getSpec());

    //return scalarMultiply(spheres[closest_sphere_index].getColor(),
                          computeLighting(P, N, -1.0 * D, spheres[closest_sphere_index].getSpec());

}

pair<double, double> InterceptRaySphere(Vector3D O, Vector3D D, Sphere sphere) {

    double r = sphere.getRadius();
    Vector3D CO = vectorSub(D,sphere.getCenter());

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

double computeLighting(Vector3D P, Vector3D N, Vector3D V, double s) {
    double i = 0.0;
    Vector3D L;
    Vector3D R;
    for (auto & light : lights) {
        if (light.getType() == "ambient") {
            i = i + light.getIntensity();
        }
        else {
            if (light.getType() == "point") {
                L = light.getVector() - P;
                //L = vectorSub(light.getVector(), P);
            } else {
                L = light.getVector();
            }
        }

        double n_dot_l = dotProduct(N,L);
        if (n_dot_l > 0.0) {
            i = i + light.getIntensity() * n_dot_l/(length(N) * length(L));
        }

        if (s != -1) {
            //R = vectorSub( 2 * dotProduct(N,L) * N, L);

            R = 2 * N * dotProduct(N, L) - L;
            //R = vectorSub(L, scalarMultiply(scalarMultiply(N, dotProduct(N,L)), 2.0));
            double r_dot_v = dotProduct(R, V);
            if (r_dot_v > 0.0) {
                i = i + light.getIntensity() * pow(r_dot_v/(length(R) * length(V)),s);
            }
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
    if (val > 255) {
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
Vector3D vectorSub(Vector3D a, Vector3D b) {
    return {a.getX()-b.getX(), a.getY()-b.getY(), a.getZ()-b.getZ()};
}


//
// Created by Adhel on 10/8/2024.
//
