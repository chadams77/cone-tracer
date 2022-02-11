
#include <GLFW/glfw3.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>

#include "vec_math.h"
#include "cl_wrapper.h"

using std::cerr;
using std::cout;
using std::endl;
using std::map;
using std::vector;
using std::unordered_map;
using std::string;
using std::stringstream;
using std::ifstream;
using std::ofstream;

const double SZ_FACTOR = 0.59;

double rayInterCount = 0.;

class VoxelNodeCL {
public:
    CLFloat4 pos_size;
    CLUByte4 normal_emit;
    CLUByte4 color_reflect;
    CLFloat offset;
    CLInt2 index_level;
    CLInt aaa, baa, aba, bba,
          aab, bab, abb, bbb;
    CLInt _dummy; // make struct size divisible by 16 bytes

    VoxelNodeCL() {
        pos_size.x = pos_size.y = pos_size.z = pos_size.w = 0.;
        normal_emit.x = normal_emit.y = normal_emit.z = normal_emit.w = 0;
        color_reflect.x = color_reflect.y = color_reflect.z = color_reflect.w = 0;
        offset = 0.;
        index_level.x = index_level.y = 0;
        aaa = baa = aba = bba = aab = bab = abb = bbb = -1;
        _dummy = 0;
    }
};

class VoxelNode {
public:
    vec3 pos, normal, color;
    double size, offset, reflect, emit;
    int level;
    int ind;

    VoxelNode * aaa = NULL;
    VoxelNode * baa = NULL;
    VoxelNode * aba = NULL;
    VoxelNode * bba = NULL;
    VoxelNode * aab = NULL;
    VoxelNode * bab = NULL;
    VoxelNode * abb = NULL;
    VoxelNode * bbb = NULL;

    VoxelNode(const vec3 & p, double s, int l=0) {
        pos = p; size = s; level=l;
        aaa = baa = aba = bba = aab = bab = abb = bbb = NULL;
        normal = vec3(0.);
        offset = 0.;
        color = vec3(0.);
        reflect = 0.;
        emit = 0.;
    }

    inline VoxelNodeCL* WriteCLNode(VoxelNodeCL * outArray, CLInt index) {
        VoxelNodeCL * out = outArray + index;
        out->pos_size.x = (CLFloat)pos.x;
        out->pos_size.y = (CLFloat)pos.y;
        out->pos_size.z = (CLFloat)pos.z;
        out->pos_size.w = (CLFloat)size;
        vec3 n2 = clamp(normal, -1., 1.);
        out->normal_emit.x = (CLUByte)floor((n2.x + 1.)*0.5*255.);
        out->normal_emit.y = (CLUByte)floor((n2.y + 1.)*0.5*255.);
        out->normal_emit.z = (CLUByte)floor((n2.z + 1.)*0.5*255.);
        out->normal_emit.w = (CLUByte)floor(clamp(emit, 0., 1.)*255.);
        vec3 c2 = clamp(color, 0., 1.);
        out->color_reflect.x = (CLUByte)floor(c2.x*255.);
        out->color_reflect.y = (CLUByte)floor(c2.y*255.);
        out->color_reflect.z = (CLUByte)floor(c2.z*255.);
        out->color_reflect.w = (CLUByte)floor(clamp(reflect, 0., 1.)*255.);
        out->offset = (CLFloat)offset;
        out->index_level.x = index;
        out->index_level.y = (CLInt)level;
        out->aaa = out->baa = out->aba = out->bba = out->aab = out->bab = out->abb = out->bbb = -1;
        out->_dummy = 0;
        ind = index;
        return out;
    }

    inline VoxelNodeCL* WriteCLNodeLinks(VoxelNodeCL * outArray, CLInt index) {
        VoxelNodeCL * out = outArray + index;
        if (aaa != NULL) { out->aaa = (CLInt)(aaa->ind); }
        if (baa != NULL) { out->baa = (CLInt)(baa->ind); }
        if (aba != NULL) { out->aba = (CLInt)(aba->ind); }
        if (bba != NULL) { out->bba = (CLInt)(bba->ind); }
        if (aab != NULL) { out->aab = (CLInt)(aab->ind); }
        if (bab != NULL) { out->bab = (CLInt)(bab->ind); }
        if (abb != NULL) { out->abb = (CLInt)(abb->ind); }
        if (bbb != NULL) { out->bbb = (CLInt)(bbb->ind); }
        return out;
    }

    inline bool rayIntersects(const vec3 & r0, const vec3 & rd, vec3 & p, double & dist, bool normalCheck = false) const {

        rayInterCount += 1.;

        vec3 bounds[2];
        bounds[0] = pos - size * 0.5;
        bounds[1] = pos + size * 0.5;

        vec3 rd2 = clamp(vec3(abs(rd.x), abs(rd.y), abs(rd.z)), 1e-9, 1e9);
        if (rd.x < 0.) {
            rd2.x = -rd2.x;
        }
        if (rd.y < 0.) {
            rd2.y = -rd2.y;
        }
        if (rd.z < 0.) {
            rd2.z = -rd2.z;
        }

        vec3 invdir = vec3(1.) / rd2;

        int sign[3];
        sign[0] = (invdir.x < 0); 
        sign[1] = (invdir.y < 0); 
        sign[2] = (invdir.z < 0); 

        double tmin, tmax, tymin, tymax, tzmin, tzmax; 
 
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

        double t;
        t = tmin;
        if (t < 0.) { 
            if (tmax < 0.) return false; 
            t = tmax;
        }
 
        dist = t;
        p = r0 + rd * dist;

        return true;
    } 

    #define _RC(nnn) if (nnn != NULL && nnn->rayIntersects(r0, rd, pTemp, distTemp) && (nodeTmp = nnn->rayCast(r0, rd, maxVoxelSize, pTemp, distTemp))) { \
            if (nodeTmp != NULL) { \
                dist = distTemp; \
                p = pTemp; \
                return nodeTmp; \
            } \
        }

    inline VoxelNode * rayCast ( const vec3 & r0, const vec3 & rd, double maxVoxelSize, vec3 & p, double & dist ) {
        if ((level == 0) || ((size/dist) <= maxVoxelSize)) {
            if (level > 0 || rayIntersects(r0, rd, p, dist, true)) {
                return this;
            }
            else {
                return NULL;
            }
        }

        VoxelNode * nodeTmp = NULL;
        vec3 pTemp = vec3(0.);
        double distTemp = 1e9;

        if (rd.z < 0.) {
            if (rd.y < 0.) {
                if (rd.x < 0.) {
                    _RC(bbb);
                    _RC(abb);
                    _RC(bab);
                    _RC(aab);
                    _RC(bba);
                    _RC(aba);
                    _RC(baa);
                    _RC(aaa);
                }
                else {
                    _RC(abb);
                    _RC(bbb);
                    _RC(aab);
                    _RC(bab);
                    _RC(aba);
                    _RC(bba);
                    _RC(aaa);
                    _RC(baa);
                }
            }
            else {
                if (rd.x < 0.) {
                    _RC(bab);
                    _RC(aab);
                    _RC(bbb);
                    _RC(abb);
                    _RC(baa);
                    _RC(aaa);
                    _RC(bba);
                    _RC(aba);
                }
                else {
                    _RC(aab);
                    _RC(bab);
                    _RC(abb);
                    _RC(bbb);
                    _RC(aaa);
                    _RC(baa);
                    _RC(aba);
                    _RC(bba);
                }
            }
        }
        else {
            if (rd.y < 0.) {
                if (rd.x < 0.) {
                    _RC(bba);
                    _RC(aba);
                    _RC(baa);
                    _RC(aaa);
                    _RC(bbb);
                    _RC(abb);
                    _RC(bab);
                    _RC(aab);
                }
                else {
                    _RC(aba);
                    _RC(bba);
                    _RC(aaa);
                    _RC(baa);
                    _RC(abb);
                    _RC(bbb);
                    _RC(aab);
                    _RC(bab);
                }
            }
            else {
                if (rd.x < 0.) {
                    _RC(baa);
                    _RC(aaa);
                    _RC(bba);
                    _RC(aba);
                    _RC(bab);
                    _RC(aab);
                    _RC(bbb);
                    _RC(abb);
                }
                else {
                    _RC(aaa);
                    _RC(baa);
                    _RC(aba);
                    _RC(bba);
                    _RC(aab);
                    _RC(bab);
                    _RC(abb);
                    _RC(bbb);
                }
            }
        }

        return NULL;
    }

};

class VResult {
public:
    VoxelNode * node;
    vec3 normal;
    vec3 color;
    vec3 p;
    double dist;

    VResult() {
        node = NULL;
    }
    VResult(const VResult & b) {
        node = b.node;
        normal = b.normal;
        color = b.color;
        p = b.p;
        dist = b.dist;
    }
};

class VLight {
public:
    vec3 pos;
    vec3 color;
    double radius;
    VLight() {
        pos = vec3(0.);
        color = vec3(0.);
        radius = 0.;
    }
    VLight(const vec3 & a, const double r, const vec3 & b) {
        pos = a;
        color = b;
        radius = r;
    }
    vec3 applyLight(const vec3 & p, const vec3 & n, const vec3 & c, const double shadowF = 1.) const {
        vec3 lightDir = pos - p;
        float dist = lightDir.length();
        lightDir /= dist;
        return (color * c) * shadowF * clamp(dot(lightDir, n), 0., 1.) * clamp(1. - (dist / radius), 0., 1.);
    }
};

class VLightCL {
public:
    CLFloat4 pos_radius;
    CLFloat4 color;

    VLightCL() {
        pos_radius.x = pos_radius.y = pos_radius.z = pos_radius.w = 0.;
        color.x = color.y = color.z = 0.;
        color.w = 1.;
    }
};

VLightCL * makeCLLightsArray(const vector<VLight> & lights, VLightCL * ret = NULL) {
    if (ret == NULL) {
        ret = new VLightCL[lights.size()];
    }
    for (size_t i=0; i<lights.size(); i++) {
        ret[i].pos_radius.x = lights[i].pos.x;
        ret[i].pos_radius.y = lights[i].pos.y;
        ret[i].pos_radius.z = lights[i].pos.z;
        ret[i].pos_radius.w = lights[i].radius;
        ret[i].color.x = lights[i].color.x;
        ret[i].color.y = lights[i].color.y;
        ret[i].color.z = lights[i].color.z;
    }
    return ret;
};

class VoxelTree {
public:
    VoxelNode * root = NULL;
    VoxelNode ** all = NULL;
    size_t nodeCount = 0;

    VoxelTree(VoxelNode * r, vector<VoxelNode*> & a) {
        nodeCount = a.size();
        all = new VoxelNode*[nodeCount];
        for (size_t i=0; i<nodeCount; i++) {
            all[i] = a[i];
        }
        root = r;
    }

    ~VoxelTree() {
        for (size_t i=0; i<nodeCount; i++) {
            delete all[i];
        }
        delete [] all;
        root = NULL;
        all = NULL;
        nodeCount = 0;
    }

    VoxelNodeCL * getVoxelsCL(VoxelNodeCL * ret = NULL) {
        if (!nodeCount) {
            return NULL;
        }
        if (ret == NULL) {
            ret = new VoxelNodeCL[nodeCount];
        }
        for (size_t index=0; index<nodeCount; index++) {
            all[index]->WriteCLNode(ret, (CLInt)index);
        }
        for (size_t index=0; index<nodeCount; index++) {
            all[index]->WriteCLNodeLinks(ret, (CLInt)index);
        }
        return ret;
    }

    inline VoxelNode * rayCast(const vec3 & r0, const vec3 & rd, double maxVoxelSize, vec3 & p, double & dist) {
        if (root != NULL && root->rayIntersects(r0, rd, p, dist)) {
            p = r0;
            dist = 0.;
            VoxelNode * ret = root->rayCast(r0, rd, maxVoxelSize, p, dist);
            if (ret == NULL) {
                p = r0;
                dist = 1e9;
            }
            else {
                p = r0 + rd * dist;
            }
            return ret;
        }
        else {
            dist = 1e9;
            p = r0 + rd * dist;
            return NULL;
        }
    }

    inline VResult rayCastCamera(
        const vec3 & cameraPosition,
        const vec3 & cameraDirection,
        const vec3 & cameraUp,
        const vec2 & screenPosition,
        const vec2 & screenSize,
        const double pixelSize,
        const vector<VLight> & lights
    ) {
        vec3 left = cross(cameraDirection, cameraUp);
        vec3 up = cross(left, cameraDirection);
        vec2 uv = (screenPosition / screenSize) - 0.5;
        vec3 r0 = left * -uv.x + up * uv.y + cameraPosition;
        vec3 r1 = (left * -uv.x + up * uv.y) * 2. + cameraPosition + cameraDirection;
        vec3 rd = normalize(r1 - r0);
        double retDist = 0.;
        vec3 retP = r0;
        VoxelNode * ret = rayCast(r0, rd, pixelSize, retP, retDist);

        VResult VR;
        VR.node = ret;
        VR.p = retP;
        VR.color = vec3(0.);
        VR.dist = retDist;
        if (ret != NULL) {
            for (size_t i=0; i<lights.size(); i++) {
                double retDist2 = 0.;
                vec3 shadowDir = normalize(lights[i].pos - retP);
                vec3 retP2 = VR.p + ret->normal * sqrt(ret->size*ret->size*8.);
                double shadowF = rayCast(vec3(retP2), shadowDir, pixelSize, retP2, retDist2) != NULL ? 0.5 : 1.;
                VR.color += lights[i].applyLight(
                    VR.p,
                    ret->normal,
                    ret->color,// * (1. - ret->reflect),
                    shadowF
                );
                if (ret->reflect > 0.) {
                    vec3 reflectDir = normalize(reflect(rd, ret->normal));
                    vec3 retP3 = retP + ret->normal * sqrt(ret->size*ret->size*8.);
                    double retDist3 = 0.;
                    VoxelNode * ret2 = rayCast(vec3(retP3), reflectDir, pixelSize, retP3, retDist3);
                    if (ret2 != NULL) {
                        double retDist4 = 0.;
                        vec3 retP4 = retP3 + ret2->normal * sqrt(ret2->size*ret2->size*8.);
                        vec3 shadowDir2 = normalize(lights[i].pos - retP4);
                        double shadowF2 = rayCast(vec3(retP4), shadowDir2, pixelSize, retP4, retDist4) != NULL ? 0.5 : 1.;
                        VR.color += lights[i].applyLight(
                            retP4,
                            ret2->normal,
                            ret2->color * ret->reflect,
                            shadowF2
                        );
                    }
                }
            }
        }
        else {
            VR.dist = 1000000.;
        }
        return VR;
    }
};

double sphereSdf(vec3 p) {
    vec3 np = normalize(p);
    return p.length() - (3000. + 1500. * pow(fnoise(np*2.5) * 0.5 + fnoise(np*5.) * 0.25 + fnoise(np*10.) * 0.125 + fnoise(np*20.) * 0.125 * 0.5, 2.));
}

VoxelTree * generateVoxelsFromSDF ( CLContext & context, string sdfName, vec3 lower, vec3 upper, double voxelSize ) {

    vec3 counts = ceil((upper - lower) / voxelSize);
    upper = lower + counts * voxelSize;

    CLInt2 windowSize;
    windowSize.x = 512;
    windowSize.y = 512;
    CLInt n = windowSize.x * windowSize.y;

    CLProgram program(context, sdfName, true);
    CLBuffer bufferOutDist(&program, n, sizeof(CLFloat), MEMORY_WRITE);
    CLBuffer bufferOutNormal(&program, n, sizeof(CLFloat3), MEMORY_WRITE);
    CLBuffer bufferOutPos(&program, n, sizeof(CLFloat3), MEMORY_WRITE);
    CLBuffer bufferOutColor(&program, n, sizeof(CLFloat3), MEMORY_WRITE);
    CLBuffer bufferOutEmit(&program, n, sizeof(CLFloat), MEMORY_WRITE);
    CLBuffer bufferOutReflect(&program, n, sizeof(CLFloat), MEMORY_WRITE);

    CLFloat3 minRange = lower, maxRange = upper;

    program.setArg("sdf_main", 0, &bufferOutDist);
    program.setArg("sdf_main", 1, &bufferOutNormal);
    program.setArg("sdf_main", 2, &bufferOutPos);
    program.setArg("sdf_main", 3, &bufferOutColor);
    program.setArg("sdf_main", 4, &bufferOutEmit);
    program.setArg("sdf_main", 5, &bufferOutReflect);
    program.setArg("sdf_main", 6, windowSize);
    program.setArg("sdf_main", 8, (CLFloat)voxelSize);
    program.setArg("sdf_main", 9, minRange);
    program.setArg("sdf_main", 10, maxRange);

    vec3 i(0., 0., 0.);
    int count = 0;
    unordered_map<int64_t, VoxelNode*> level[2];
    level[0] = unordered_map<int64_t, VoxelNode*>();
    level[1] = unordered_map<int64_t, VoxelNode*>();
    vector<VoxelNode *> allVN;
    key3 key = key3(0, 0, 0);
    double SF = sqrt(voxelSize*voxelSize*4.);
    for (i.z=lower.z+voxelSize*0.5, key.z=0; i.z <= (upper.z+0.1); i.z += voxelSize, key.z++) {
        for (i.x=lower.x+voxelSize*0.5, key.x=0; i.x <= (upper.x+0.1); i.x += voxelSize * (double)windowSize.x, key.x += windowSize.x) {
            for (i.y=lower.y+voxelSize*0.5, key.y=0; i.y <= (upper.y+0.1); i.y += voxelSize * (double)windowSize.y, key.y += windowSize.y) {

                program.setArg("sdf_main", 7, CLInt3(key.x, key.y, key.z));
                if (!program.callFunction("sdf_main", n)) {
                    exit(0);
                }

                bufferOutDist.readSync();
                bufferOutNormal.readSync();
                bufferOutPos.readSync();
                bufferOutColor.readSync();
                bufferOutEmit.readSync();
                bufferOutReflect.readSync();

                CLFloat * aDist = bufferOutDist.dataFloat();
                CLFloat3 * aPos = bufferOutPos.dataFloat3();
                CLFloat3 * aNormal = bufferOutNormal.dataFloat3();
                CLFloat3 * aClr = bufferOutColor.dataFloat3();
                CLFloat * aEmit = bufferOutEmit.dataFloat();
                CLFloat * aReflect = bufferOutReflect.dataFloat();

                for (int ox=0; ox<windowSize.x; ox++) {
                    for (int oy=0; oy<windowSize.y; oy++) {
                        int off = ox + oy * windowSize.x;
                        CLFloat dist = aDist[off];
                        CLFloat3 pos = aPos[off];
                        if (dist > -SF*2. && pos.x <= (upper.x+0.1) && pos.y <= (upper.y+0.1) && pos.z <= (upper.z+0.1)) {
                            count += 1;
                            vec3 p(pos.x, pos.y, pos.z);
                            CLFloat3 normal = aNormal[off];
                            VoxelNode * vn = new VoxelNode(p, voxelSize, 0);
                            vn->color = vec3(aClr[off].x, aClr[off].y, aClr[off].z);
                            vn->reflect = aReflect[off];
                            vn->emit = aEmit[off];
                            vn->normal = vec3(normal.x, normal.y, normal.z);
                            vn->offset = dist;
                            key3 k2(key.x + ox, key.y + oy, key.z);
                            level[0][k2.hash()] = vn;
                            allVN.push_back(vn);
                        }
                    }
                }
            }
        }
        cerr << "completed slice: " << key.z << ", count: " << count << endl;
    }
    cerr << "voxel count: " << count << endl;

    int last = 0, next = 1;
    while (level[last].size() > 1) {
        cerr << "generating level " << (level[last].begin()->second->level + 1) << "... ";
        level[next].clear();
        key3 minKey(100000, 100000, 100000);
        key3 maxKey(-100000, -100000, -100000);
        for (unordered_map<int64_t, VoxelNode*>::iterator ii=level[last].begin(); ii!=level[last].end(); ii++) {
            key3 k = key3(ii->first);
            minKey.x = k.x < minKey.x ? k.x : minKey.x;
            minKey.y = k.y < minKey.y ? k.y : minKey.y;
            minKey.z = k.z < minKey.z ? k.z : minKey.z;
            maxKey.x = k.x > maxKey.x ? k.x : maxKey.x;
            maxKey.y = k.y > maxKey.y ? k.y : maxKey.y;
            maxKey.z = k.z > maxKey.z ? k.z : maxKey.z;
        }
        minKey = minKey / 2;
        maxKey = maxKey / 2;
        for (int x=minKey.x; x<=maxKey.x; x++) {
            for (int y=minKey.y; y<=maxKey.y; y++) {
                for (int z=minKey.z; z<=maxKey.z; z++) {
                    key3 k2 = key3(x, y, z) * 2;
                    double chCount = 0.;
                    VoxelNode * parent = new VoxelNode(vec3(0.), 0., -1);
                    parent->color = vec3(0.);
                    parent->emit = 0.;
                    parent->reflect = 0.;
                    for (int xo=0; xo<=1; xo++) {
                        for (int yo=0; yo<=1; yo++) {
                            for (int zo=0; zo<=1; zo++) {
                                key3 k3 = k2 + key3(xo, yo, zo);
                                unordered_map<int64_t, VoxelNode*>::iterator fi = level[last].find(k3.hash());
                                if (fi != level[last].end()) {
                                    VoxelNode * child = fi->second;
                                    chCount += 1.;
                                    if (parent->level < 0)  {
                                        parent->level = child->level + 1;
                                        parent->size = child->size * 2.;
                                    }
                                    parent->normal += child->normal;
                                    parent->offset += child->offset;
                                    parent->color += child->color;
                                    parent->emit += child->emit;
                                    parent->reflect += child->reflect;
                                    if (xo < 0.5) { // a
                                        if (yo < 0.5) { // a
                                            if (zo < 0.5) { //a
                                                parent->aaa = child;
                                                parent->pos = child->pos + vec3(1., 1., 1.) * child->size * 0.5;
                                            }
                                            else { // b
                                                parent->aab = child;
                                                parent->pos = child->pos + vec3(1., 1., -1.) * child->size * 0.5;
                                            }
                                        }
                                        else { // b
                                            if (zo < 0.5) { //a
                                                parent->aba = child;
                                                parent->pos = child->pos + vec3(1., -1., 1.) * child->size * 0.5;
                                            }
                                            else { // b
                                                parent->abb = child;
                                                parent->pos = child->pos + vec3(1., -1., -1.) * child->size * 0.5;
                                            }
                                        }
                                    }
                                    else { // b
                                        if (yo < 0.5) { // a
                                            if (zo < 0.5) { //a
                                                parent->baa = child;
                                                parent->pos = child->pos + vec3(-1., 1., 1.) * child->size * 0.5;
                                            }
                                            else { // b
                                                parent->bab = child;
                                                parent->pos = child->pos + vec3(-1., 1., -1.) * child->size * 0.5;
                                            }
                                        }
                                        else { // b
                                            if (zo < 0.5) { //a
                                                parent->bba = child;
                                                parent->pos = child->pos + vec3(-1., -1., 1.) * child->size * 0.5;
                                            }
                                            else { // b
                                                parent->bbb = child;
                                                parent->pos = child->pos + vec3(-1., -1., -1.) * child->size * 0.5;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (chCount > 0.0001) {
                        key3 key2 = k2 / 2;
                        parent->normal = normalize(parent->normal);
                        parent->color /= chCount;
                        parent->emit /= chCount;
                        parent->reflect /= chCount;
                        level[next][key2.hash()] = parent;
                        allVN.push_back(parent);
                    }
                    else {
                        delete parent;
                        parent = NULL;
                    }
                }
            }
        }
        last = 1 - last;
        next = 1 - next;
        cerr << "node count: " << level[last].size() << endl;
    }

    VoxelTree * ret = new VoxelTree(level[last].begin()->second, allVN);

    return ret;

}

int main(void)
{
    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit()) {
        return -1;
    }

    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(SCREEN_HEIGHT, SCREEN_HEIGHT, "Cone Tracing Test", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    CLContext cl_context;

    CLFloat voxelSize = 10000./512.;

    vector<VLight> lights;
    lights.push_back(VLight(vec3(-10000., -10000., -10000.), 30000., vec3(2., 2., 2.)));
    lights.push_back(VLight(vec3(0., 0., 0.), 2000., vec3(2., 0., 0.)));

    VoxelTree * tree = generateVoxelsFromSDF(cl_context, "box_inception", vec3(-5000.), vec3(5000.), voxelSize);

    CLProgram program(cl_context, "raytrace", false);
    CLBuffer vTreeBfr(&program, tree->nodeCount, sizeof(VoxelNodeCL), MEMORY_READ);
    CLBuffer vLightBfr(&program, lights.size(), sizeof(VLightCL), MEMORY_READ);

    CLImageGL outImage(&program, SCREEN_WIDTH, SCREEN_HEIGHT, MEMORY_WRITE);
    CLBuffer depthBfr(&program, SCREEN_WIDTH * SCREEN_HEIGHT, sizeof(CLFloat), MEMORY_WRITE);

    tree->getVoxelsCL((VoxelNodeCL*)(vTreeBfr.data));
    makeCLLightsArray(lights, (VLightCL*)(vLightBfr.data));

    vTreeBfr.writeSync();
    vLightBfr.writeSync();

    program.setArg("raytrace_main", 0, &vTreeBfr);
    program.setArg("raytrace_main", 1, &vLightBfr);
    program.setArg("raytrace_main", 2, &outImage);
    program.setArg("raytrace_main", 3, &depthBfr);

    CLInt rootIndex = tree->root->ind;

    program.setArg("raytrace_main", 8, rootIndex);
    program.setArg("raytrace_main", 9, (CLFloat)voxelSize);
    program.setArg("raytrace_main", 10, (CLFloat)20000.f);
    program.setArg("raytrace_main", 11, (CLInt)lights.size());

    //depthBfr.readSync();

    int index = 0;

    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window)) {

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0.0, 1., 0.0, 1., -1.0, 1.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        double ang = ((double)index/60) * 3.141592 * 0.25;
        double r = (sin(ang*1.5f) * 0.5 + 0.5) * 8000. + 2000.;
        vec3 camPos = vec3(cos(ang) * r, 0., sin(ang) * r);
        vec3 camDir = normalize(vec3(0.) - camPos);
        vec3 camUp = normalize(vec3(0., 1., 0.));

        CLInt2 screenSize(SCREEN_WIDTH, SCREEN_HEIGHT);

        program.setArg("raytrace_main", 4, screenSize);
        program.setArg("raytrace_main", 5, CLFloat3(camDir));
        program.setArg("raytrace_main", 6, CLFloat3(camPos));
        program.setArg("raytrace_main", 7, CLFloat3(camUp));

        program.acquireImageGL(&outImage);

        if (!program.callFunction("raytrace_main", SCREEN_WIDTH * SCREEN_HEIGHT)) {
            exit(0);
        }

        program.releaseImageGL(&outImage);

        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, outImage.glTex);

        glBegin(GL_QUADS);
            glTexCoord2f(0., 0.); glVertex2f(0., 0.);
            glTexCoord2f(0., 1.); glVertex2f(0., 1.);
            glTexCoord2f(1., 1.); glVertex2f(1., 1.);
            glTexCoord2f(1., 0.); glVertex2f(1., 0.);
        glEnd();

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();

        index += 1;
        if (!(index % 60)) {
            cerr << "60 frames." << endl;
        }
    }

    delete tree;
    tree = NULL;

    glfwTerminate();
    return 0;
}