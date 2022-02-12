typedef struct __attribute__((packed)) _VoxelNode {
    float4 pos_size;
    uchar4 normal_emit;
    uchar4 color_reflect;
    float offset;
    int2 index_level;
    int aaa, baa, aba, bba,
        aab, bab, abb, bbb;
    int _dummy; // make struct size divisible by 16 bytes
} VoxelNode;

typedef struct __attribute__((packed)) _LightSource {
    float4 pos_range;
    float4 color;
} LightSource;

__global VoxelNode *VOX;
__global LightSource *LIGHT;
float3 gCamPos;
float gVoxSize;
int gLightCount;

bool rayIntersects(__global VoxelNode * node, float3 r0, float3 rd, float3 * p, float * dist, float gVMul) {

    float3 bounds[2];

    gVMul = 0.f;

    float3 sz = (float3)((node->pos_size.w + gVoxSize * gVMul) * 0.5f);

    bounds[0] = node->pos_size.xyz - sz;
    bounds[1] = node->pos_size.xyz + sz;

    float3 rd2 = clamp((float3)(fabs(rd.x), fabs(rd.y), fabs(rd.z)), 1e-9f, 1e9f);
    if (rd.x < 0.f) {
        rd2.x = -rd2.x;
    }
    if (rd.y < 0.f) {
        rd2.y = -rd2.y;
    }
    if (rd.z < 0.f) {
        rd2.z = -rd2.z;
    }

    float3 invdir = (float3)(1.f) / rd2;

    int sign[3];
    sign[0] = (invdir.x < 0.f); 
    sign[1] = (invdir.y < 0.f); 
    sign[2] = (invdir.z < 0.f); 

    float tmin, tmax, tymin, tymax, tzmin, tzmax; 
 
    tmin = (bounds[sign[0]].x - r0.x) * invdir.x; 
    tmax = (bounds[1-sign[0]].x - r0.x) * invdir.x; 
    tymin = (bounds[sign[1]].y - r0.y) * invdir.y; 
    tymax = (bounds[1-sign[1]].y - r0.y) * invdir.y; 
 
    if ((tmin > tymax) || (tymin > tmax)) {
        return false; 
    }
    if (tymin > tmin) {
        tmin = tymin; 
    }
    if (tymax < tmax) {
        tmax = tymax; 
    }
 
    tzmin = (bounds[sign[2]].z - r0.z) * invdir.z; 
    tzmax = (bounds[1-sign[2]].z - r0.z) * invdir.z; 

    if ((tmin > tzmax) || (tzmin > tmax)) {
        return false; 
    }
    if (tzmin > tmin) {
        tmin = tzmin; 
    }

    if (tzmax < tmax) {
        tmax = tzmax; 
    }

    float t;
    t = tmin;
    if (t < 0.f) { 
        if (tmax < 0.f) return false; 
        t = tmax;
    }

    *dist = t;
    *p = r0 + rd * t;

    return true;
}

#define STACK_SIZE 20
#define MAX_ITS 256
#define _RC(nnn) if ( (nnn) >= 0 && n >= stack[stackPtr].y) { \
                    __global VoxelNode * nodeTmp = VOX + (nnn); \
                    if (rayIntersects(nodeTmp, r0, rd, &pTemp, &distTemp, 1.0)) { \
                        *p = pTemp; \
                        *dist = distTemp; \
                        stack[stackPtr].y = n+1; \
                        stackPtr += 1; \
                        stack[stackPtr] = (int2)(nnn, 0); \
                        continue; \
                    } \
                } \
                n += 1; \
                stack[stackPtr].y = max(stack[stackPtr].y, n);
__global VoxelNode * rayCast(int rootIndex, float3 r0, float3 rd, float maxVoxelSize, float3 * p, float * dist, float maxDist) {
    int2 stack[STACK_SIZE];
    int stackPtr = 0;
    stack[stackPtr] = (int2)(rootIndex, 0);
    if (rayIntersects(VOX + rootIndex, r0, rd, p, dist, 1.0)) {
        if ((*dist) > maxDist) {
            return NULL;
        }
        *p = r0;
        *dist = 0.f;

        float3 pTemp = r0;
        float distTemp = 0.f;

        int its = 0;

        while (stackPtr >= 0 && its < MAX_ITS) {
            __global VoxelNode * node = VOX + stack[stackPtr].x;
            if ((*dist) > maxDist) {
                return NULL;
            }
            if (node->index_level.y <= 0 || ((node->pos_size.w/distance(*p, gCamPos)) <= maxVoxelSize)) {
                return node;
            }

            int n = 0;

            if (rd.z < 0.f) {
                if (rd.y < 0.f) {
                    if (rd.x < 0.f) {
                        _RC(node->bbb)
                        _RC(node->abb)
                        _RC(node->bab)
                        _RC(node->aab)
                        _RC(node->bba)
                        _RC(node->aba)
                        _RC(node->baa)
                        _RC(node->aaa)
                    }
                    else {
                        _RC(node->abb)
                        _RC(node->bbb)
                        _RC(node->aab)
                        _RC(node->bab)
                        _RC(node->aba)
                        _RC(node->bba)
                        _RC(node->aaa)
                        _RC(node->baa)
                    }
                }
                else {
                    if (rd.x < 0.f) {
                        _RC(node->bab)
                        _RC(node->aab)
                        _RC(node->bbb)
                        _RC(node->abb)
                        _RC(node->baa)
                        _RC(node->aaa)
                        _RC(node->bba)
                        _RC(node->aba)
                    }
                    else {
                        _RC(node->aab)
                        _RC(node->bab)
                        _RC(node->abb)
                        _RC(node->bbb)
                        _RC(node->aaa)
                        _RC(node->baa)
                        _RC(node->aba)
                        _RC(node->bba)
                    }
                }
            }
            else {
                if (rd.y < 0.f) {
                    if (rd.x < 0.f) {
                        _RC(node->bba)
                        _RC(node->aba)
                        _RC(node->baa)
                        _RC(node->aaa)
                        _RC(node->bbb)
                        _RC(node->abb)
                        _RC(node->bab)
                        _RC(node->aab)
                    }
                    else {
                        _RC(node->aba)
                        _RC(node->bba)
                        _RC(node->aaa)
                        _RC(node->baa)
                        _RC(node->abb)
                        _RC(node->bbb)
                        _RC(node->aab)
                        _RC(node->bab)
                    }
                }
                else {
                    if (rd.x < 0.f) {
                        _RC(node->baa)
                        _RC(node->aaa)
                        _RC(node->bba)
                        _RC(node->aba)
                        _RC(node->bab)
                        _RC(node->aab)
                        _RC(node->bbb)
                        _RC(node->abb)
                    }
                    else {
                        _RC(node->aaa)
                        _RC(node->baa)
                        _RC(node->aba)
                        _RC(node->bba)
                        _RC(node->aab)
                        _RC(node->bab)
                        _RC(node->abb)
                        _RC(node->bbb)
                    }
                }
            }

            stackPtr -= 1;
            its += 1;
        }
    }
    else {
        *dist = 1e9f;
        *p = r0 + rd * (*dist);
    }
    return NULL;
}

#define MAX_LIGHTS 16
float4 parseNormalEmit(__global VoxelNode * n) {
    float4 ret;
    ret.x = ((float)(n->normal_emit.x)) / 255.f;
    ret.y = ((float)(n->normal_emit.y)) / 255.f;
    ret.z = ((float)(n->normal_emit.z)) / 255.f;
    ret.xyz = normalize((ret.xyz - 0.5f)*2.f);
    ret.w = ((float)(n->normal_emit.w)) / 255.f;
    return ret;
}
float4 parseColorReflect(__global VoxelNode * n) {
    float4 ret;
    ret.x = ((float)(n->color_reflect.x)) / 255.f;
    ret.y = ((float)(n->color_reflect.y)) / 255.f;
    ret.z = ((float)(n->color_reflect.z)) / 255.f;
    ret.w = ((float)(n->color_reflect.w)) / 255.f;
    return ret;
}

float3 reflect(float3 I, float3 N) {
    return I - (float3)(2.0f * dot(N, I)) * N;
}

#define NUM_REFLECT 1

float3 getColor(int rootIndex, float3 r0, float3 rd, float maxVoxelSize, float maxDist) {
    
    float3 p;
    float dist;

    float3 ret = (float3)(0., 0., 0.);
    float lastReflect = 1.f;

    for (int ri=0; ri<=NUM_REFLECT; ri++) {

        __global VoxelNode * N1 = rayCast(rootIndex, r0, rd, maxVoxelSize * (ri > 0 ? 2.f : 1.f), &p, &dist, maxDist);

        if (N1 == NULL) {
            break;
        }

        float F1 = sqrt(N1->pos_size.w*N1->pos_size.w*8.f);

        float4 ce1 = parseColorReflect(N1);
        float4 nr1 = parseNormalEmit(N1);

        ret += ce1.xyz * nr1.w * (1.f - ce1.w) * lastReflect;

        for (int i=0; i<gLightCount && i<MAX_LIGHTS; i++) {
            __global LightSource * L = LIGHT + i;
            float distF = distance(p, L->pos_range.xyz);
            if (distF < L->pos_range.w) {
                float3 lightDir = (L->pos_range.xyz - p) / distF;
                float3 _t1; float _t2;
                float shadowF = rayCast(rootIndex, p + nr1.xyz * F1, lightDir, maxVoxelSize*2.f, &_t1, &_t2, L->pos_range.w) == NULL ? 1.f : 0.25f;
                distF = 1.f - distF / L->pos_range.w;
                ret += ce1.xyz * distF * L->color.xyz * clamp(dot(lightDir, nr1.xyz), 0.f, 1.f) * shadowF * (1.f - ce1.w) * lastReflect;
            }
        }

        if (ce1.w <= 0.f) {
            break;
        }

        r0 = p + nr1.xyz * F1;
        rd = normalize(reflect(rd, nr1.xyz));

        lastReflect = ce1.w;
    }

    return ret;

}

__kernel void raytrace_main( __global VoxelNode *voxels,
                             __global LightSource *lights,
                             __write_only image2d_t out_color,
                             __global float  *out_depth,
                             int2 size,
                             float3 cameraDirection,
                             float3 cameraPosition,
                             float3 cameraUp,
                             int rootIndex,
                             float voxelSize,
                             float maxDist,
                             int lightCount)
{
    int id = get_global_id(0);
    int n = size.x * size.y;

    if (id < n) {

        VOX = voxels;
        LIGHT = lights;
        gLightCount = lightCount;

        int x = id % size.x;
        int y = (id - x) / size.x;

        float2 screenPosition = (float2)((float)x, (float)y);
        float2 screenSize = (float2)((float)size.x, (float)size.y);

        gCamPos = cameraPosition;
        gVoxSize = voxelSize;
        float maxVoxelSize = 0.2f * voxelSize / max(screenSize.x, screenSize.y);

        float3 left = cross(cameraDirection, cameraUp);
        float3 up = cross(left, cameraDirection);
        float2 uv = (screenPosition / screenSize) - 0.5f;
        float3 r0 = left * -uv.x + up * uv.y + cameraPosition;
        float3 r1 = (left * -uv.x + up * uv.y) * 2.f + cameraPosition + cameraDirection;
        float3 rd = normalize(r1 - r0);

        float4 clr = (float4)(0., 0., 0., 1.);

        clr.xyz = getColor(rootIndex, r0, rd, maxVoxelSize, maxDist);

        clr = clamp(clr, (float4)(0.), (float4)(1.));
        write_imagef(out_color, (int2)(x, y), clr);

    }
}