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

class object {
public:
    object(const Vector3D &P0, const Vector3D &P1, const Vector3D &P2, const Vector3D &color, const double spec, const double reflect) {
        this->name = "Triangle";
        this->P0 = P0;
        this->P1 = P1;
        this->P2 = P2;
        this->color = color;
        this->spec = spec;
        this->reflect = reflect;
        this->radius = 0;
    }
    object(const Vector3D &center, const double radius, const Vector3D &color, const double spec, const double reflect) {
        this->name = "Sphere";
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
    [[nodiscard]] Vector3D getP0() const {
        return P0;
    }
    [[nodiscard]] Vector3D getP1() const {
        return P1;
    }
    [[nodiscard]] Vector3D getP2() const {
        return P2;
    }
    [[nodiscard]] string getName() const {
        return name;
    }

private:
    string name;
    Vector3D center;
    double radius;
    Vector3D color;
    double spec;
    double reflect;
    Vector3D P0;
    Vector3D P1;
    Vector3D P2;
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
constexpr double inf = 9999999.0;
constexpr double d = 1;
constexpr int recursion_depth = 3;
Vector3D origin(0.0, 0.0, 0.0);
Vector3D bColor(0.0, 0.0, 0.0);

vector<object> objects;
vector<Light> lights;

void imageGen();
Vector3D canvasToViewport(double x,double y);
Vector3D traceRay(const Vector3D &O, const Vector3D &D, double t_min, double t_max, int recursion_depth);
pair<double, double> InterceptRaySphere(const Vector3D &O, const Vector3D &D, const object &sphere);
pair<double, double> InterceptRayTriangle(const Vector3D &O, const Vector3D &D, const object &Triangle);
double computeLighting(const Vector3D& P, const Vector3D& N, const Vector3D& V, double s);
bool inRange(const double &min, const double &max, const double &value);
double length(const Vector3D &vec);
void addValues();
double dotProduct(const Vector3D &a, const Vector3D &b);
Vector3D crossProduct(const Vector3D &a, const Vector3D &b);
int max(double val);
pair<int, double> closestIntersection(const Vector3D &O, const Vector3D &D, double t_min, double t_max);
Vector3D ReflectRay(const Vector3D &R, const Vector3D &N);
bool pointInTriangle(const Vector3D& P, Vector3D A, Vector3D B, Vector3D C);

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

    const int closest_index = fst;
    const double closest_t = snd;

    if (closest_index == -1) {
        return bColor;
    }
    const Vector3D P = O + closest_t * D;
    Vector3D N = P - objects[closest_index].getCenter();
    N = N/length(N);

    const Vector3D localColor = objects[closest_index].getColor() * computeLighting(P, N, -1*D, objects[closest_index].getSpec());

    double r = objects[closest_index].getReflect();
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
    pair<double, double> values;
    for (int i = 0; i < objects.size(); i++) {
        if (objects[i].getName() == "Sphere") {
            values = InterceptRaySphere(O,D,objects[i]);
        } else {
            values = InterceptRayTriangle(O,D,objects[i]);
        }

        const double t1 = values.first;
        const double t2 = values.second;

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

pair<double, double> InterceptRaySphere(const Vector3D &O, const Vector3D &D, const object &sphere) {

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

pair<double, double> InterceptRayTriangle(const Vector3D &O, const Vector3D &D, const object &Triangle) {

    Vector3D P0 = Triangle.getP0();
    Vector3D P1 = Triangle.getP1();
    Vector3D P2 = Triangle.getP2();

    Vector3D edge1 = P1 - P0;
    Vector3D edge2 = P2 - P0;
    Vector3D normal = crossProduct(edge1, edge2);

    double NdotD = dotProduct(normal, D);
    if (NdotD == 0) {
        return {inf, inf};
    }

    double t = dotProduct(normal, P0 - O) / NdotD;
    if (t < 0) {
        return {inf, inf};
    }

    Vector3D P = O + D * t;

    if (pointInTriangle(P, P0, P1, P2)) {
        return {t, inf};
    } else {
        return {inf, inf};
    }
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

bool pointInTriangle(const Vector3D &P, Vector3D A, Vector3D B, Vector3D C) {
    Vector3D v0 = C - A;
    Vector3D v1 = B - A;
    Vector3D v2 = P - A;

    double dot00 = dotProduct(v0, v0);
    double dot01 = dotProduct(v0, v1);
    double dot02 = dotProduct(v0, v2);
    double dot11 = dotProduct(v1, v1);
    double dot12 = dotProduct(v1, v2);

    double invDenom = 1.0/(dot00*dot11-dot01*dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    return (u >= 0) && (v >= 0) && (u + v <= 1);
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
Vector3D crossProduct(const Vector3D &a, const Vector3D &b) {
    double x = a.getY() * b.getZ() - a.getZ() * b.getY();
    double y = a.getZ() * b.getX() - a.getX() * b.getZ();
    double z = a.getX() * b.getY() - a.getY() * b.getX();

    return {x, y, z};
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
    const auto red = object(Vector3D(0,-1,3), 1.0, Vector3D(255,0,0), 500, 0.2);
    const auto blue = object(Vector3D(-2,0,4), 1.0, Vector3D(0,0,255), 10, 0.3);
    const auto green = object(Vector3D(2,0,4), 1.0, Vector3D(0,255,0), 500, 0.4);
    // const auto yellow = object(Vector3D(0,-5001,0), 5000, Vector3D(255,255,0), 1000, 0.5);

    const auto triangle = object(Vector3D(0,-1,3), Vector3D(-2,0,4), Vector3D(2,0,4), Vector3D(255,255,0), 500, 0.5);

    objects.push_back(red);
    objects.push_back(blue);
    objects.push_back(green);
    //objects.push_back(yellow);
    objects.push_back(triangle);

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
