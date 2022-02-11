#pragma once

class vec3 {
public:
    double x, y, z;
    inline vec3() { x=y=z=0.; }
    inline vec3(double _x, double _y, double _z) { x=_x; y=_y; z=_z; }
    inline vec3(double _x) { x=_x; y=_x; z=_x; }
    inline vec3(const vec3 & b) { x=b.x; y=b.y; z=b.z; }
    inline vec3 operator+(const vec3 & b) const { return vec3(x+b.x, y+b.y, z+b.z); }
    inline vec3 operator-(const vec3 & b) const { return vec3(x-b.x, y-b.y, z-b.z); }
    inline vec3 operator+(const double & b) const { return vec3(x+b, y+b, z+b); }
    inline vec3 operator-(const double & b) const { return vec3(x-b, y-b, z-b); }
    inline vec3 operator*(const vec3 & b) const { return vec3(x*b.x, y*b.y, z*b.z); }
    inline vec3 operator/(const vec3 & b) const { return vec3(x/b.x, y/b.y, z/b.z); }
    inline vec3 operator*(const double & b) const { return vec3(x*b, y*b, z*b); }
    inline vec3 operator/(const double & b) const { return vec3(x/b, y/b, z/b); }
    inline vec3 & operator+=(const vec3 & b) { x += b.x; y += b.y; z += b.z; return *this; }
    inline vec3 & operator-=(const vec3 & b) { x -= b.x; y -= b.y; z -= b.z; return *this; }
    inline vec3 & operator+=(const double & b) { x += b; y += b; z += b; return *this; }
    inline vec3 & operator-=(const double & b) { x -= b; y -= b; z -= b; return *this; }
    inline vec3 & operator*=(const vec3 & b) { x *= b.x; y *= b.y; z *= b.z; return *this; }
    inline vec3 & operator/=(const vec3 & b) { x /= b.x; y /= b.y; z /= b.z; return *this; }
    inline vec3 & operator*=(const double & b) { x *= b; y *= b; z *= b; return *this; }
    inline vec3 & operator/=(const double & b) { x /= b; y /= b; z /= b; return *this; }
    inline double length() const { return sqrt(x*x+y*y+z*z); }
    inline double distance(const vec3 & b) const { return (b - (*this)).length(); }
    inline bool operator==(const vec3 & b) const {
        return abs(x-b.x) < 1e-9 && abs(y-b.y) < 1e-9 && abs(z-b.z) < 1e-9;
    }
    inline bool operator<(const vec3 & b) const { 
        vec3 delta = (*this) - b;
        if (delta.x < -1e-9) {
            return true;
        }
        else if (delta.y > 1e-9) {
            return false;
        }
        else {
            if (delta.y < -1e-9) {
                return true;
            }
            else if (delta.y > 1e-9) {
                return false;
            }
            else {
                return delta.z < -1e-9;
            }
        }
    }
};

class vec2 {
public:
    double x, y;
    inline vec2() { x=y=0.; }
    inline vec2(double _x=0., double _y=0.) { x=_x; y=_y; }
    inline vec2(double _x) { x=_x; y=_x; }
    inline vec2(const vec2 & b) { x=b.x; y=b.y; }
    inline vec2 operator+(const vec2 & b) const { return vec2(x+b.x, y+b.y); }
    inline vec2 operator-(const vec2 & b) const { return vec2(x-b.x, y-b.y); }
    inline vec2 operator+(const double & b) const { return vec2(x+b, y+b); }
    inline vec2 operator-(const double & b) const { return vec2(x-b, y-b); }
    inline vec2 operator*(const vec2 & b) const { return vec2(x*b.x, y*b.y); }
    inline vec2 operator/(const vec2 & b) const { return vec2(x/b.x, y/b.y); }
    inline vec2 operator*(const double & b) const { return vec2(x*b, y*b); }
    inline vec2 operator/(const double & b) const { return vec2(x/b, y/b); }
    inline vec2 & operator+=(const vec2 & b) { x += b.x; y += b.y; return *this; }
    inline vec2 & operator-=(const vec2 & b) { x -= b.x; y -= b.y; return *this; }
    inline vec2 & operator+=(const double & b) { x += b; y += b; return *this; }
    inline vec2 & operator-=(const double & b) { x -= b; y -= b; return *this; }
    inline vec2 & operator*=(const vec2 & b) { x *= b.x; y *= b.y; return *this; }
    inline vec2 & operator/=(const vec2 & b) { x /= b.x; y /= b.y; return *this; }
    inline vec2 & operator*=(const double & b) { x *= b; y *= b; return *this; }
    inline vec2 & operator/=(const double & b) { x /= b; y /= b; return *this; }
    inline double length() const { return sqrt(x*x+y*y); }
    inline double distance(const vec2 & b) const { return (b - (*this)).length(); }
    inline bool operator<(const vec2 & b) const { 
        if (x < b.x) {
            return true;
        }
        else if (x > b.x) {
            return false;
        }
        else {
            return y < b.y;
        }
    }
};

class key3 {
public:
    int x, y, z;
    inline key3() { x=y=z= 0; }
    inline key3(int a, int b, int c) { x=a; y=b; z=c; }
    inline key3(const key3 & b) { x=b.x; y=b.y; z=b.z; }
    inline key3(int64_t a) {
        z = (a & ((1LL<<20LL)-1LL)) - (1LL << 19LL);
        y = ((a >> 20LL) & ((1LL<<20LL)-1LL)) - (1LL << 19LL);
        x = (a >> 40LL) - (1LL << 19LL);
    }
    inline key3 operator+(const key3 & b) const { return key3(x+b.x, y+b.y, z+b.z); }
    inline key3 operator-(const key3 & b) const { return key3(x-b.x, y-b.y, z-b.z); }
    inline key3 operator+(const int & b) const { return key3(x+b, y+b, z+b); }
    inline key3 operator-(const int & b) const { return key3(x-b, y-b, z-b); }
    inline key3 operator*(const key3 & b) const { return key3(x*b.x, y*b.y, z*b.z); }
    inline key3 operator/(const key3 & b) const { return key3(x/b.x, y/b.y, z/b.z); }
    inline key3 operator*(const int & b) const { return key3(x*b, y*b, z*b); }
    inline key3 operator/(const int & b) const { return key3(x/b, y/b, z/b); }
    inline bool operator==(const key3 & b) const {
        return (x == b.x && y == b.y && z == b.z);
    }
    inline bool operator<(const key3 & b) const {
        if (x == b.x) {
            if (y == b.y) {
                return z < b.z;
            }
            return y < b.y;
        }
        return x < b.x;
    }
    inline int64_t hash() const {
        return (((int64_t)x + (1LL << 19LL)) << 40LL) + (((int64_t)y + (1LL << 19LL)) << 20LL) + ((int64_t)z + (1LL << 19LL));
    }
};

static inline vec3 clamp (const vec3 & v, const vec3 & a, const vec3 & b) {
    return vec3(
        v.x < a.x ? a.x : v.x > b.x ? b.x : v.x,
        v.y < a.y ? a.y : v.y > b.y ? b.y : v.y,
        v.z < a.z ? a.z : v.z > b.z ? b.z : v.z
    );
}

static inline vec2 clamp (const vec2 & v, const vec2 & a, const vec2 & b) {
    return vec2(
        v.x < a.x ? a.x : v.x > b.x ? b.x : v.x,
        v.y < a.y ? a.y : v.y > b.y ? b.y : v.y
    );
}

static inline vec3 clamp (const vec3 & v, const double a, const double b) {
    return vec3(
        v.x < a ? a : v.x > b ? b : v.x,
        v.y < a ? a : v.y > b ? b : v.y,
        v.z < a ? a : v.z > b ? b : v.z
    );
}

static inline vec2 clamp (const vec2 & v, const double a, const double b) {
    return vec2(
        v.x < a ? a : v.x > b ? b : v.x,
        v.y < a ? a : v.y > b ? b : v.y
    );
}

static inline double clamp (double v, const double a, const double b) {
    return v < a ? a : v > b ? b : v;
}

static inline vec3 floor (const vec3 & v) {
    return vec3(
        floor(v.x),
        floor(v.y),
        floor(v.z)
    );
}

static inline vec2 floor (const vec2 & v) {
    return vec2(
        floor(v.x),
        floor(v.y)
    );
}

static inline vec3 fract (const vec3 & v) {
    return v - floor(v);
}

static inline vec2 fract (const vec2 & v) {
    return v - floor(v);
}

static inline double fract (const double v) {
    return v - floor(v);
}

static inline vec3 ceil (const vec3 & v) {
    return vec3(
        ceil(v.x),
        ceil(v.y),
        ceil(v.z)
    );
}

static inline vec2 ceil (const vec2 & v) {
    return vec2(
        ceil(v.x),
        ceil(v.y)
    );
}

static inline vec3 round (const vec3 & v) {
    return vec3(
        round(v.x),
        round(v.y),
        round(v.z)
    );
}

static inline vec2 round (const vec2 & v) {
    return vec2(
        round(v.x),
        round(v.y)
    );
}

static inline vec3 mix (const vec3 & a, const vec3 & b, const double t) {
    if (t < 0.) {
        return a;
    }
    else if (t > 1.) {
        return b;
    }
    else {
        return a + (b - a) * t;
    }
}

static inline vec2 mix (const vec2 & a, const vec2 & b, const double t) {
    if (t < 0.) {
        return a;
    }
    else if (t > 1.) {
        return b;
    }
    else {
        return a + (b - a) * t;
    }
}

static inline double mix (const double a, const double b, const double t) {
    if (t < 0.) {
        return a;
    }
    else if (t > 1.) {
        return b;
    }
    else {
        return a + (b - a) * t;
    }
}

static inline double dot (const vec3 & a, const vec3 & b) {
    vec3 t = a * b;
    return t.x + t.y + t.z;
}

static inline double dot (const vec2 & a, const vec2 & b) {
    vec2 t = a * b;
    return t.x + t.y;
}

static inline vec3 cross (const vec3 & a, const vec3 & b) {
    return vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

static inline vec3 reflect (const vec3 & I, const vec3 & N) {
    return I - vec3(2.0 * dot(N, I)) * N;
}

static inline double cross (const vec2 & a, const vec2 & b) {
    return a.x * b.y - a.y * b.x;
}

static inline vec3 normalize (const vec3 & v) {
    return v / v.length();
}

static inline vec2 normalize (const vec2 & v) {
    return v / v.length();
}

static inline float fhash(vec3 p) {
    p = fract( p*0.3183099+.1 );
    p *= 17.0;
    return fract( p.x*p.y*p.z*(p.x+p.y+p.z) );
}

static inline float fnoise(vec3 x) {
    vec3 i = floor(x);
    vec3 f = fract(x);
    f = f*f*(f*-2.0 + 3.0);
    
    return mix(mix(mix( fhash(i+vec3(0,0,0)), 
                        fhash(i+vec3(1,0,0)),f.x),
                    mix( fhash(i+vec3(0,1,0)), 
                        fhash(i+vec3(1,1,0)),f.x),f.y),
                mix(mix( fhash(i+vec3(0,0,1)), 
                        fhash(i+vec3(1,0,1)),f.x),
                    mix( fhash(i+vec3(0,1,1)), 
                        fhash(i+vec3(1,1,1)),f.x),f.y),f.z);
}