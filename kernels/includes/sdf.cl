__kernel void sdf_main( __global float  *out_dist,
                        __global float3 *out_normal,                   
                        __global float3 *out_pos,
                        __global float3 *out_color,
                        __global float *out_emissive,
                        __global float *out_reflective,
                        int2 size,
                        int3 offset,
                        float voxelSize,
                        float3 minRange,
                        float3 maxRange)
{
    int id = get_global_id(0);
    int n = size.x * size.y;

    if (id < n) {

        int x = id % size.x;
        int y = (id - x) / size.x;
        int z = offset.z;
        x += offset.x;
        y += offset.y;
        float3 p = (float3)(x, y, z);
        p = p * (float3)(voxelSize);
        p += minRange + (float3)(voxelSize*0.5);

        float SF = sqrt(voxelSize*voxelSize*4.);
        out_color[id] = (float3)(0.f, 0.f, 0.f);
        out_emissive[id] = out_reflective[id] = 0.;
        float dist = sdf(p, out_color + id, out_emissive + id, out_reflective + id);
        if (dist < SF) {
            float nDist = voxelSize;
            float3 _a; float _b; float _c;
            float A = sdf(p + (float3)(1.,0.,0.) * nDist, &_a, &_b, &_c), B = sdf(p - (float3)(1.,0.,0.) * nDist, &_a, &_b, &_c),
                  C = sdf(p + (float3)(0.,1.,0.) * nDist, &_a, &_b, &_c), D = sdf(p - (float3)(0.,1.,0.) * nDist, &_a, &_b, &_c),
                  E = sdf(p + (float3)(0.,0.,1.) * nDist, &_a, &_b, &_c), F = sdf(p - (float3)(0.,0.,1.) * nDist, &_a, &_b, &_c);
            float3 normal = (float3)(
                A - B,
                C - D,
                E - F
            );
            float nlength = length(normal);
            if (nlength < 0.00001) {
                normal = (float3)(0., 0., -1.);
            }
            else {
                normal /= nlength;
            }
            if (sdf(p + normal * SF * 2.f, &_a, &_b, &_c) < 0.f) {
                normal = -normal;
            }
            if (A > SF || B > SF || C > SF || D > SF || E > SF || F > SF) {
                out_normal[id] = normal;
                out_dist[id] = max(dist, -SF*0.9f);
            }
            else {
                out_normal[id] = (float3)(0., 0., 0.);
                out_dist[id] = -100000.;
            }
        }
        else {
            out_normal[id] = (float3)(0., 0., 0.);
            out_dist[id] = -100000.;
        }

        out_pos[id] = p;
    }
}