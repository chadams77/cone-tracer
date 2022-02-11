float sdf(const float3 _p, float3 * color, float * emissive, float * reflective) {
    float3 p = _p;//rotate3D(_p, (float3)(1., 1., 1.), (float)M_PI * 0.5f);
    float dist = 100000.f;
    *color = (float3)(0.5, 0.5, 0.5);
    *emissive = 0.;
    *reflective = 0.5;
    for (int level=1; level<=5; level++) {
        float size = 4000.f / pow(2.f, (float)level);
        float3 p2 = floor(p / size);
        int3 ip2 = (int3)(p2.x, p2.y, p2.z);
        if (fabs(p2.x) > pow(2.f, (float)level - 1.f) || fabs(p2.y) > pow(2.f, (float)level - 1.f) || fabs(p2.z) > pow(2.f, (float)level - 1.f)) {
            continue;
        }
        if ((ip2.x % 2) && (ip2.y % 3) && !(ip2.z % 4)) {
            float3 c = (p2 + 0.5f) * size;
            //dist = min(dist, length(p - c) - 0.4f * 0.95f * size);
            dist = min(dist, sdBoxSigned(p, c, 0.5f * 0.95f * (float3)size) - fnoise(normalize(p-c)*size*0.02f) * size * 0.025f);
            break;
        }
    }
    return dist;
}