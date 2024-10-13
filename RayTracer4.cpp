#include <iostream>
#include <cmath>
#include <fstream>
#include <optional>
#include <utility>
#include <vector>

using namespace std;

class Vector3D {
public:

    Vector3D(const double &x, const double &y, const double &z) {
        v[0] = x;
        v[1] = y;
        v[2] = z;
    }

    Vector3D() : v{0,0,0} { }

    [[nodiscard]] double getX() const {
        return v[0];
    }

    [[nodiscard]] double getY() const {
        return v[1];
    }

    [[nodiscard]] double getZ() const {
        return v[2];
    }

    friend ostream& operator <<(ostream& out, const Vector3D& obj){
        out << "<" << obj.v[0] << "," << obj.v[1] << "," << obj.v[2] << ">";

        return out;
    }

    friend Vector3D operator *(const double scalar, const Vector3D& obj){
        Vector3D tmp;

        tmp.v[0] = scalar*obj.v[0];
        tmp.v[1] = scalar*obj.v[1];
        tmp.v[2] = scalar*obj.v[2];

        return tmp;
    }
    friend Vector3D operator *(const Vector3D& obj, const double scalar){
        Vector3D tmp;

        tmp.v[0] = obj.v[0]*scalar;
        tmp.v[1] = obj.v[1]*scalar;
        tmp.v[2] = obj.v[2]*scalar;

        return tmp;
    }
    friend Vector3D operator *(const Vector3D& a, const Vector3D& b) {
        Vector3D tmp;

        tmp.v[0] = a.v[0]*b.v[0];
        tmp.v[1] = a.v[1]*b.v[1];
        tmp.v[2] = a.v[2]*b.v[2];

        return tmp;
    }

    friend Vector3D operator /(const double scalar, const Vector3D& obj) {
        Vector3D tmp;

        tmp.v[0] = scalar/obj.v[0];
        tmp.v[1] = scalar/obj.v[1];
        tmp.v[2] = scalar/obj.v[2];

        return tmp;
    }
    friend Vector3D operator /(const Vector3D& obj,const double scalar) {
        Vector3D tmp;

        tmp.v[0] = obj.v[0]/scalar;
        tmp.v[1] = obj.v[1]/scalar;
        tmp.v[2] = obj.v[2]/scalar;

        return tmp;
    }
    friend Vector3D operator /(const Vector3D& a, const Vector3D& b) {
        Vector3D tmp;

        tmp.v[0] = a.v[0]/b.v[0];
        tmp.v[1] = a.v[1]/b.v[1];
        tmp.v[2] = a.v[2]/b.v[2];

        return tmp;
    }

    friend Vector3D operator -(const double& scalar, const Vector3D& obj) {
        Vector3D tmp;

        tmp.v[0] = obj.v[0]-scalar;
        tmp.v[1] = obj.v[1]-scalar;
        tmp.v[2] = obj.v[2]-scalar;

        return tmp;
    }
    friend Vector3D operator -(const Vector3D& obj, const double& scalar) {
        Vector3D tmp;

        tmp.v[0] = scalar-obj.v[0];
        tmp.v[1] = scalar-obj.v[1];
        tmp.v[2] = scalar-obj.v[2];

        return tmp;
    }
    friend Vector3D operator -(const Vector3D& a, const Vector3D& b) {
        Vector3D tmp;

        tmp.v[0] = a.v[0]-b.v[0];
        tmp.v[1] = a.v[1]-b.v[1];
        tmp.v[2] = a.v[2]-b.v[2];

        return tmp;
    }

    friend Vector3D operator +(const double& scalar, const Vector3D& obj) {
        Vector3D tmp;

        tmp.v[0] = obj.v[0]+scalar;
        tmp.v[1] = obj.v[1]+scalar;
        tmp.v[2] = obj.v[2]+scalar;

        return tmp;
    }
    friend Vector3D operator +(const Vector3D& obj, const double& scalar) {
        Vector3D tmp;

        tmp.v[0] = scalar+obj.v[0];
        tmp.v[1] = scalar+obj.v[1];
        tmp.v[2] = scalar+obj.v[2];

        return tmp;
    }
    friend Vector3D operator +(const Vector3D& a, const Vector3D& b) {
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
    Sphere(const Vector3D &center, const double radius, const Vector3D &color, const double spec, const double reflect) {
        this->center = center;
        this->radius = radius;
        this->color = color;
        this->spec = spec;
        this->reflect = reflect;
    }

    [[nodiscard]] Vector3D getCenter() const {
        return center;
    }

    [[nodiscard]] double getRadius() const {
        return radius;
    }

    [[nodiscard]] Vector3D getColor() const {
        return color;
    }
    [[nodiscard]] double getSpec() const {
        return spec;
    }

    [[nodiscard]] double getReflect() const {
        return reflect;
    }

private:
    Vector3D center;
    double radius;
    Vector3D color;
    double spec;
    double reflect;
};

class Light {
public:
    Light(const double intensity, string type, const Vector3D &vec) {
        this->intensity = intensity;
        this->type = std::move(type);
        this->vec = vec;
    }

    [[nodiscard]] double getIntensity() const {return intensity;}
    [[nodiscard]] Vector3D getVector() const {return vec;}
    string getType() {return type;}

private:
    double intensity;
    Vector3D vec;
    string type;
};

constexpr int Cw = 1000;
constexpr int Ch = 1000;
constexpr int Vw = 1;
constexpr int Vh = 1;
constexpr int inf = 9999999;
constexpr double d = 1;
constexpr int recursion_depth = 3;
Vector3D origin(0.0, 0.0, 0.0);
Vector3D bColor(0.0, 0.0, 0.0);

vector<Sphere> spheres;
vector<Light> lights;

void imageGen();
Vector3D canvasToViewport(double x,double y);
Vector3D traceRay(const Vector3D &O, const Vector3D &D, double t_min, double t_max, int recursion_depth);
pair<double, double> InterceptRaySphere(const Vector3D &O, const Vector3D &D, const Sphere &sphere);
double computeLighting(const Vector3D& P, const Vector3D& N, const Vector3D& V, double s);
bool inRange(const double &min, const double &max, const double &value);
double length(const Vector3D &vec);
void addValues();
double dotProduct(const Vector3D &a, const Vector3D &b);
int max(double val);
pair<int, double> closestIntersection(const Vector3D &O, const Vector3D &D, double t_min, double t_max);
Vector3D ReflectRay(const Vector3D &R, const Vector3D &N);

Vector3D vectorAdd(const Vector3D &a, const Vector3D &b);
Vector3D vectorSub(const Vector3D &a, const Vector3D &b);

Vector3D scalarMultiply(const Vector3D &a, double b);
Vector3D scalarDivide(const Vector3D &a, double b);


int main(int argc, char const *argv[]) {
    addValues();
    imageGen();
    return 0;
}
Vector3D canvasToViewport(const double x, const double y) {
    return {x*Vw/Cw, y*Vh/Ch, d};
}

Vector3D traceRay(const Vector3D &O, const Vector3D &D, const double t_min, const double t_max, const int recursion_depth) {

    auto [fst, snd] = closestIntersection(O, D, t_min, t_max);

    const int closest_sphere_index = fst;
    const double closest_t = snd;

    if (closest_sphere_index == -1) {
        return bColor;
    }
    const Vector3D P = O + closest_t * D;
    Vector3D N = P - spheres[closest_sphere_index].getCenter();
    N = N/length(N);

    const Vector3D localColor = spheres[closest_sphere_index].getColor() * computeLighting(P, N, -1*D, spheres[closest_sphere_index].getSpec());

    double r = spheres[closest_sphere_index].getReflect();
    if (recursion_depth <= 0 || r <= 0) {
        return localColor;
    }

    const Vector3D R = ReflectRay(-1*D, N);
    const Vector3D reflected_color = traceRay(P, R, 0.001, t_max, recursion_depth - 1);

    return localColor * (1 - r) + reflected_color * r;
}

pair<int, double> closestIntersection(const Vector3D &O, const Vector3D &D, const double t_min, const double t_max) {

    double closest_t = inf;
    int closest_sphere_index = -1;
    for (int i = 0; i < spheres.size(); i++) {
        auto [fst, snd] = InterceptRaySphere(O,D,spheres[i]);
        const double t1 = fst;
        const double t2 = snd;

        if (inRange(t_min, t_max, t1) && t1 < closest_t) {
            closest_t = t1;
            closest_sphere_index = i;
        }

        if (inRange(t_min, t_max, t2) && t2 < closest_t) {
            closest_t = t2;
            closest_sphere_index = i;
        }
    }

    return {closest_sphere_index, closest_t};
}

pair<double, double> InterceptRaySphere(const Vector3D &O, const Vector3D &D, const Sphere &sphere) {

    const double r = sphere.getRadius();
    const Vector3D CO = vectorSub(O,sphere.getCenter());

    const double a = dotProduct(D,D);
    const double b = 2 * dotProduct(D,CO);
    const double c = dotProduct(CO,CO) - (double)(r*r);

    const double dis = b*b - 4*a*c;

    if (dis < 0) {
        return pair<double, double>{inf,inf};
    }

    double t1 = (-b + sqrt(dis)) / (2*a);
    double t2 = (-b - sqrt(dis)) / (2*a);

    return { t1,t2 };
}

double computeLighting(const Vector3D& P, const Vector3D& N, const Vector3D& V, const double s) {
    double i = 0.0;
    double t_max = 0.0;
    Vector3D L;

    for (auto & light : lights) {
        if (light.getType() == "ambient") {
            i += light.getIntensity();
        }
        else {
            if (light.getType() == "point") {
                L = light.getVector() - P;
                t_max = 1.0;
            } else {
                L = light.getVector();
                t_max = 100000.0;
            }
        }

        auto [fst, snd] = closestIntersection(P, L, 0.001, t_max);
        if (const int index = fst; index != -1) {
            continue;
        }

        if (const double n_dot_l = dotProduct(N, L); n_dot_l > 0.0) {
            i += light.getIntensity() * n_dot_l / (length(N) * length(L));
        }

        if (s != -1) {
            const Vector3D R = 2 * N * dotProduct(N, L) - L;
            if (const double r_dot_v = dotProduct(R, V); r_dot_v > 0.0) {
                i += light.getIntensity() * pow(r_dot_v / (length(R) * length(V)), s);
            }
        }
    }
    return i;
}

Vector3D ReflectRay(const Vector3D &R, const Vector3D &N) {
    return 2 * N * dotProduct(N, R) - R;
}








double length(const Vector3D &vec) {
    return sqrt(vec.getX()*vec.getX()+vec.getY()*vec.getY()+vec.getZ()*vec.getZ());
}

bool inRange(const double &min, const double &max, const double &value) {
    if (value > min && value < max) {
        return true;
    }
    return false;
}

int max(const double val) {
    if (val > 255) {
        return 255;
    }
    if (val < 0) {
        return 0;
    }
    return static_cast<int>(val);
}


Vector3D scalarMultiply(const Vector3D &a, const double b) {
    return {a.getX() * b, a.getY() * b, a.getZ() * b};
}
Vector3D scalarDivide(const Vector3D &a, const double b) {
    return {a.getX() / b, a.getY() / b, a.getZ() / b};
}


Vector3D vectorAdd(const Vector3D &a, const Vector3D &b) {
    return {a.getX()+b.getX(), a.getY()+b.getY(), a.getZ()+b.getZ()};
}
Vector3D vectorSub(const Vector3D &a, const Vector3D &b) {
    return {a.getX()-b.getX(), a.getY()-b.getY(), a.getZ()-b.getZ()};
}

double dotProduct(const Vector3D &a, const Vector3D &b) {
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
            Vector3D color = traceRay(origin, D, d, inf, recursion_depth);

            int r = max(color.getX());
            int g = max(color.getY());
            int b = max(color.getZ());

            ppmFile << r << " " << g << " " << b << " ";
        }
        ppmFile << "\n";
    }

    ppmFile.close();

    cout << "PPM file created: output.ppm" << endl;
}

void addValues() {
    const auto red = Sphere(Vector3D(0,-1,3), 1.0, Vector3D(255,0,0), 500, 0.2);
    const auto blue = Sphere(Vector3D(-2,0,4), 1.0, Vector3D(0,0,255), 10, 0.3);
    const auto green = Sphere(Vector3D(2,0,4), 1.0, Vector3D(0,255,0), 500, 0.4);
    const auto yellow = Sphere(Vector3D(0,-5001,0), 5000, Vector3D(255,255,0), 1000, 0.5);
    spheres.push_back(red);
    spheres.push_back(blue);
    spheres.push_back(green);
    spheres.push_back(yellow);

    const auto ambient = Light(0.2, "ambient",Vector3D(0,0,0));
    const auto point = Light(0.6, "point", Vector3D(2,1,0));
    const auto directional = Light(0.2, "directional", Vector3D(1,4,4));
    lights.push_back(ambient);
    lights.push_back(point);
    lights.push_back(directional);

}

//
// Created by Adhel on 10/8/2024.
//
