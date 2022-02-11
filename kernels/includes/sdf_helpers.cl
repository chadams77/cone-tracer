float fhash(const float3 _p) {
    float3 _dummy;
    float3 p = fract(_p*0.3183099f+.1f, &_dummy);
    p *= 17.0f;
    float _dummy2;
    return fract(p.x*p.y*p.z*(p.x+p.y+p.z), &_dummy2);
}

float fnoise(float3 x) {
    float3 _dummy;
    float3 i = floor(x);
    float3 f = fract(x, &_dummy);
    f = f*f*(f*-2.0f + 3.0f);
    
    return mix(mix(mix( fhash(i+(float3)(0.f,0.f,0.f)), 
                        fhash(i+(float3)(1.f,0.f,0.f)),f.x),
                    mix( fhash(i+(float3)(0.f,1.f,0.f)), 
                        fhash(i+(float3)(1.f,1.f,0.f)),f.x),f.y),
                mix(mix( fhash(i+(float3)(0.f,0.f,1.f)), 
                        fhash(i+(float3)(1.f,0.f,1.f)),f.x),
                    mix( fhash(i+(float3)(0.f,1.f,1.f)), 
                        fhash(i+(float3)(1.f,1.f,1.f)),f.x),f.y),f.z);
}

float sdBoxSigned (const float3 p, const float3 c, const float3 s) {
    float x = max(
        p.x - c.x - s.x * 0.5,
        c.x - p.x - s.x * 0.5
    );

    float y = max(
        p.y - c.y - s.y * 0.5,
        c.y - p.y - s.y * 0.5
    );

    float z = max(
        p.z - c.z - s.z * 0.5,
        c.z - p.z - s.z * 0.5
    );

    float d = x;
    d = max(d, y);
    d = max(d, z);
    return d;
}

float3 rotate3D(const float3 p, const float3 _axis, const float angle) {

    float3 axis = normalize(_axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;

    float4 p4 = (float4)(p.x, p.y, p.z, 1.);

    float4 m4_1 = (float4)(oc * axis.x * axis.x + c, oc * axis.x * axis.y - axis.z * s, oc * axis.z * axis.x + axis.y * s, 0.0);
    float4 m4_2 = (float4)(oc * axis.x * axis.y + axis.z * s, oc * axis.y * axis.y + c, oc * axis.y * axis.z - axis.x * s, 0.0);
    float4 m4_3 = (float4)(oc * axis.z * axis.x - axis.y * s, oc * axis.y * axis.z + axis.x * s, oc * axis.z * axis.z + c, 0.0);
    float4 m4_4 = (float4)(0.0, 0.0, 0.0, 1.0);

    float4 p4_2 = (float4)(
        p4.x*m4_1.x + p4.y*m4_1.y + p4.z*m4_1.z + p4.w*m4_1.w,
        p4.x*m4_2.x + p4.y*m4_2.y + p4.z*m4_2.z + p4.w*m4_2.w,
        p4.x*m4_3.x + p4.y*m4_3.y + p4.z*m4_3.z + p4.w*m4_3.w,
        p4.x*m4_4.x + p4.y*m4_4.y + p4.z*m4_4.z + p4.w*m4_4.w
    );

    return p4_2.xyz;

}