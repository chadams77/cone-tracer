float sdf(const float3 p, float3 * color, float * emissive, float * reflective) {
    float3 np = normalize(p);
    float lev = pow(fnoise(np*2.5f) * 0.5f + fnoise(np*5.f) * 0.25f + fnoise(np*10.f) * 0.125f + fnoise(np*20.f) * 0.125f * 0.5f + fnoise(np*40.f) * 0.125f * 0.25f + fnoise(np*80.f) * 0.125f * 0.125f, 2.f);

    *color = mix((float3)(0.2f, 0.2f, 0.2f), (float3)(0.4f, 1.f, 0.4f), lev);
    *emissive = pow(lev, 4.f);
    *reflective = 0.1f;

    if (lev < 0.25f) {
        lev = 0.25f;
        *reflective = 0.9f;
    }

    return length(p) - (3000.f + 1500.f * lev);
}