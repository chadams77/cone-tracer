#pragma once

#include <GLFW/glfw3.h>
#include <CL/cl.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>

using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::map;
using std::stringstream;
using std::ifstream;
using std::ofstream;

#include "vec_math.h"

const GLenum PIXEL_FORMAT = GL_RGBA;
GLuint SCREEN_WIDTH = 1024;
GLuint SCREEN_HEIGHT = 1024;
GLuint DATA_SIZE = SCREEN_WIDTH * SCREEN_HEIGHT * 4;
GLubyte * colorBuffer = NULL;

enum DeviceType 
{ 
	DEVICE_CPU			= (1 << 1), 
	DEVICE_GPU			= (1 << 2)
};

enum MemoryType
{
	MEMORY_READ			= (1 << 2),
	MEMORY_WRITE		= (1 << 1),
	MEMORY_READ_WRITE	= (1 << 0)
};

#define CLFloat float_t
#define CLDouble double_t
#define CLByte int8_t
#define CLShort int16_t
#define CLInt int32_t
#define CLLong int64_t
#define CLUByte uint8_t
#define CLUShort uint16_t
#define CLUInt uint32_t
#define CLULong uint64_t

#pragma pack(push)

#pragma pack(1)
class CLFloat2 {
public:
    union { CLFloat x, r; };
    union { CLFloat y, g; };
    CLFloat2() {
        x = y = 0.;
    }
    CLFloat2(const vec2 & v) {
        x = v.x;
        y = v.y;
    }
};

#pragma pack(1)
class CLFloat3 {
public:
    union { CLFloat x, r; };
    union { CLFloat y, g; };
    union { CLFloat z, b; };
    CLFloat _dummy;
    CLFloat3() {
        x = y = z = 0.;
    }
    CLFloat3(const vec3 & v) {
        x = v.x;
        y = v.y;
        z = v.z;
    }
};

#pragma pack(1)
class CLFloat4 {
public:
    union { CLFloat x, r; };
    union { CLFloat y, g; };
    union { CLFloat z, b; };
    union { CLFloat w, a; };
    CLFloat4() {
        x = y = z = w = 0.;
    }
    CLFloat4(const vec3 & v, double _w = 0.) {
        x = v.x;
        y = v.y;
        z = v.z;
        w = _w;
    }
};

#pragma pack(1)
class CLInt2 {
public:
    union { CLInt x, r; };
    union { CLInt y, g; };
    CLInt2() { x = y = 0; }
    CLInt2(CLInt _x, CLInt _y) { x = _x; y = _y; }
};

#pragma pack(1)
class CLInt3 {
public:
    union { CLInt x, r; };
    union { CLInt y, g; };
    union { CLInt z, b; };
    CLInt3() { x = y = z = 0; }
    CLInt3(CLInt _x, CLInt _y, CLInt _z) { x = _x; y = _y; z = _z; }
    CLInt _dummy;
};

#pragma pack(1)
class CLInt4 {
public:
    union { CLInt x, r; };
    union { CLInt y, g; };
    union { CLInt z, b; };
    union { CLInt w, a; };
};

#pragma pack(1)
class CLUInt2 {
public:
    union { CLUInt x, r; };
    union { CLUInt y, g; };
};

#pragma pack(1)
class CLUInt3 {
public:
    union { CLUInt x, r; };
    union { CLUInt y, g; };
    union { CLUInt z, b; };
    CLUInt _dummy;
};

#pragma pack(1)
class CLUInt4 {
public:
    union { CLUInt x, r; };
    union { CLUInt y, g; };
    union { CLUInt z, b; };
    union { CLUInt w, a; };
};

#pragma pack(1)
class CLDouble2 {
public:
    union { CLDouble x, r; };
    union { CLDouble y, g; };
};

#pragma pack(1)
class CLDouble3 {
public:
    union { CLDouble x, r; };
    union { CLDouble y, g; };
    union { CLDouble z, b; };
    CLDouble _dummy;
};

#pragma pack(1)
class CLDouble4 {
public:
    union { CLDouble x, r; };
    union { CLDouble y, g; };
    union { CLDouble z, b; };
    union { CLDouble w, a; };
};

#pragma pack(1)
class CLUByte4 {
public:
    union { CLUByte x, r; };
    union { CLUByte y, g; };
    union { CLUByte z, b; };
    union { CLUByte w, a; };
};

#pragma pack(pop)

class CLContext {
public:
    cl::Context				context;
    cl::CommandQueue		queue;
    vector<cl::Platform>	platforms;
    vector<cl::Device>		devices;
    size_t                  preferredDevice;
    size_t                  preferredPlatform;
    size_t                  preferredDeviceWorkload;

    CLContext() {

        cl::Platform::get(&platforms);

        preferredDevice = 0;
        preferredPlatform = 0;
        preferredDeviceWorkload = 0;

        for (size_t j=0; j<platforms.size(); j++) {
            platforms[j].getDevices(static_cast<cl_device_type>(DEVICE_GPU), &devices);
            for (size_t i=0; i<devices.size(); i++) {
                size_t workload = devices[i].getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
                if (workload > preferredDeviceWorkload) {
                    preferredDevice = i;
                    preferredPlatform = j;
                    preferredDeviceWorkload = workload;
                }
            }
        }

        platforms[preferredPlatform].getDevices(static_cast<cl_device_type>(DEVICE_GPU), &devices);

        cerr << "OpenGL/CL Context\n" << "Name: " << devices[preferredDevice].getInfo<CL_DEVICE_NAME>()
            << "\nVendor: " << devices[preferredDevice].getInfo<CL_DEVICE_VENDOR>() 
            << "\nDriver Version: " << devices[preferredDevice].getInfo<CL_DRIVER_VERSION>() 
            << "\nDevice Profile: " << devices[preferredDevice].getInfo<CL_DEVICE_PROFILE>() 
            << "\nDevice Version: " << devices[preferredDevice].getInfo<CL_DEVICE_VERSION>()
            << "\nMax Work Group Size: " << devices[preferredDevice].getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>()
            << endl;
        
        cl_context_properties properties[] = {
            CL_GL_CONTEXT_KHR, (cl_context_properties)wglGetCurrentContext(),
            CL_WGL_HDC_KHR, (cl_context_properties)wglGetCurrentDC(),
            CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[preferredPlatform])(),
            0
        };

        context = cl::Context(devices, properties);

    }

    void ReportError(CLInt err, string prefix) {
        switch (err) {
            case CL_SUCCESS:
                return;
                break;
            case CL_DEVICE_NOT_FOUND:
                cerr << prefix << "Device not found" << endl;
                break;
            case CL_DEVICE_NOT_AVAILABLE:
                cerr << prefix << "Device not available" << endl;
                break;
            case CL_COMPILER_NOT_AVAILABLE:
                cerr << prefix << "Compiler not available" << endl;
                break;
            case CL_MEM_OBJECT_ALLOCATION_FAILURE:
                cerr << prefix << "Memory object allocation failure" << endl;
                break;
            case CL_OUT_OF_RESOURCES:
                cerr << prefix << "Out of resources" << endl;
                break;
            case CL_OUT_OF_HOST_MEMORY:
                cerr << prefix << "Out of host memory" << endl;
                break;
            case CL_PROFILING_INFO_NOT_AVAILABLE:
                cerr << prefix << "Profiling info not available" << endl;
                break;
            case CL_MEM_COPY_OVERLAP:
                cerr << prefix << "Memory copy overlap" << endl;
                break;
            case CL_IMAGE_FORMAT_MISMATCH:
                cerr << prefix << "Image format mismatch" << endl;
                break;
            case CL_IMAGE_FORMAT_NOT_SUPPORTED:
                cerr << prefix << "Image format not supported" << endl;
                break;
            case CL_BUILD_PROGRAM_FAILURE:
                cerr << prefix << "Build program failure" << endl;
                break;
            case CL_MAP_FAILURE:
                cerr << prefix << "Map failure" << endl;
                break;
            case CL_INVALID_VALUE:
                cerr << prefix << "Invalid value" << endl;
                break;
            case CL_INVALID_DEVICE_TYPE:
                cerr << prefix << "Invalid device type" << endl;
                break;
            case CL_INVALID_PLATFORM:
                cerr << prefix << "Invalid platform" << endl;
                break;
            case CL_INVALID_DEVICE:
                cerr << prefix << "Invalid device" << endl;
                break;
            case CL_INVALID_CONTEXT:
                cerr << prefix << "Invalid context" << endl;
                break;
            case CL_INVALID_QUEUE_PROPERTIES:
                cerr << prefix << "Invalid queue properties" << endl;
                break;
            case CL_INVALID_COMMAND_QUEUE:
                cerr << prefix << "Invalid command queue" << endl;
                break;
            case CL_INVALID_HOST_PTR:
                cerr << prefix << "Invalid host pointer" << endl;
                break;
            case CL_INVALID_MEM_OBJECT:
                cerr << prefix << "Invalid memory object" << endl;
                break;
            case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
                cerr << prefix << "Invalid image format descriptor" << endl;
                break;
            case CL_INVALID_IMAGE_SIZE:
                cerr << prefix << "Invalid image size" << endl;
                break;
            case CL_INVALID_SAMPLER:
                cerr << prefix << "Invalid sampler" << endl;
                break;
            case CL_INVALID_BINARY:
                cerr << prefix << "Invalid binary" << endl;
                break;
            case CL_INVALID_BUILD_OPTIONS:
                cerr << prefix << "Invalid build options" << endl;
                break;
            case CL_INVALID_PROGRAM:
                cerr << prefix << "Invalid program" << endl;
                break;
            case CL_INVALID_PROGRAM_EXECUTABLE:
                cerr << prefix << "Invalid program executable" << endl;
                break;
            case CL_INVALID_KERNEL_NAME:
                cerr << prefix << "Invalid kernel name" << endl;
                break;
            case CL_INVALID_KERNEL_DEFINITION:
                cerr << prefix << "Invalid kernel definition" << endl;
                break;
            case CL_INVALID_KERNEL:
                cerr << prefix << "Invalid kernel" << endl;
                break;
            case CL_INVALID_ARG_INDEX:
                cerr << prefix << "Invalid argument index" << endl;
                break;
            case CL_INVALID_ARG_VALUE:
                cerr << prefix << "Invalid argument value" << endl;
                break;
            case CL_INVALID_ARG_SIZE:
                cerr << prefix << "Invalid argument size" << endl;
                break;
            case CL_INVALID_KERNEL_ARGS:
                cerr << prefix << "Invalid kernel arguments" << endl;
                break;
            case CL_INVALID_WORK_DIMENSION:
                cerr << prefix << "Invalid work dimension" << endl;
                break;
            case CL_INVALID_WORK_GROUP_SIZE:
                cerr << prefix << "Invalid work group size" << endl;
                break;
            case CL_INVALID_WORK_ITEM_SIZE:
                cerr << prefix << "invalid work item size" << endl;
                break;
            case CL_INVALID_GLOBAL_OFFSET:
                cerr << prefix << "Invalid global offset" << endl;
                break;
            case CL_INVALID_EVENT_WAIT_LIST:
                cerr << prefix << "Invalid event wait list" << endl;
                break;
            case CL_INVALID_EVENT:
                cerr << prefix << "Invalid event" << endl;
                break;
            case CL_INVALID_OPERATION:
                cerr << prefix << "Invalid operation" << endl;
                break;
            case CL_INVALID_GL_OBJECT:
                cerr << prefix << "Invalid OpenGL object" << endl;
                break;
            case CL_INVALID_BUFFER_SIZE:
                cerr << prefix << "Invalid buffer size" << endl;
                break;
            case CL_INVALID_MIP_LEVEL:
                cerr << prefix << "Invalid MIP level" << endl;
                break;
            default:
                break;
        }
        return;
    }
};

class CLBuffer;
class CLBufferGL;
class CLImageGL;

class CLProgram {
public:
    cl::Program::Sources source;
    cl::Program program;
    cl_int err = CL_SUCCESS;
    cl::Event event;
    cl::CommandQueue queue;
    map<string, cl::Kernel*> functions;
    CLContext * context;

    CLProgram(CLContext & _context, string filename, bool isSDF = false) {
        context = &_context;
        ifstream file((string("kernels/") + filename + (isSDF ? "_sdf.cl" : ".cl")).c_str());
        stringstream buffer;
        buffer << file.rdbuf();
        string code = buffer.str();

        if (isSDF) {
            ifstream file2("kernels/includes/sdf_helpers.cl");
            stringstream buffer2;
            buffer2 << file2.rdbuf();
            code = buffer2.str() + string("\n") + code;
        }

        if (isSDF) {
            ifstream file2("kernels/includes/sdf.cl");
            stringstream buffer2;
            buffer2 << file2.rdbuf();
            code = code + string("\n") + buffer2.str();
        }

        err = CL_SUCCESS;
        source = cl::Program::Sources(1, std::make_pair(code.c_str(), code.length()));
        context->ReportError(err, filename + ": ");
        program = cl::Program(context->context, source);
        program.build(context->devices);
        queue = cl::CommandQueue(context->context, context->devices[context->preferredDevice], 0, &err);
        context->ReportError(err, filename + ": ");

        cl::STRING_CLASS errStr = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(context->devices[context->preferredDevice]);
        if (errStr.length() > 1 && errStr.find("error") != string::npos) {
            cerr << filename << ":\n" << errStr << endl;
            exit(0);
        }
    }

    cl::Kernel * getFunction(string function) {
        if (functions.find(function) != functions.end()) {
            return functions[function];
        }
        cl::Kernel * kernel = new cl::Kernel(program, function.c_str(), &err);
        context->ReportError(err, "getFunction: ");
        functions[function] = kernel;
        return kernel;
    }

    bool callFunction(string function, size_t n) {
        cl::Kernel * kernel = getFunction(function);
        size_t mwSize = kernel->getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(context->devices[context->preferredDevice]);
        size_t mul = kernel->getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(context->devices[context->preferredDevice]);
        size_t localSize = ((size_t)floor((float)mwSize / (float)mul)) * mul;
        size_t globalSize = ((n+localSize-1) / localSize) * localSize;
        cl::Event event;
        int err = queue.enqueueNDRangeKernel(*kernel, cl::NullRange, cl::NDRange(globalSize), cl::NDRange(localSize), NULL, &event);
        context->ReportError(err, "callFunction: ");
        event.wait();
        return err == CL_SUCCESS;
    }

    #define _SET_ARG(_TYPE) void setArg(string function, int arg, _TYPE v) { \
        cl_int err = getFunction(function)->setArg(arg, sizeof(_TYPE), &v); \
        if (err != CL_SUCCESS) { \
            stringstream ss; \
            ss << "(" << arg << "): "; \
            context->ReportError(err, ss.str()); \
        } \
    }

    _SET_ARG(CLFloat)
    _SET_ARG(CLDouble)
    _SET_ARG(CLByte)
    _SET_ARG(CLShort)
    _SET_ARG(CLInt)
    _SET_ARG(CLLong)
    _SET_ARG(CLUByte)
    _SET_ARG(CLUShort)
    _SET_ARG(CLUInt)
    _SET_ARG(CLULong)
    _SET_ARG(CLFloat2)
    _SET_ARG(CLFloat3)
    _SET_ARG(CLFloat4)
    _SET_ARG(CLInt2)
    _SET_ARG(CLInt3)
    _SET_ARG(CLInt4)
    _SET_ARG(CLUInt2)
    _SET_ARG(CLUInt3)
    _SET_ARG(CLUInt4)
    _SET_ARG(CLDouble2)
    _SET_ARG(CLDouble3)
    _SET_ARG(CLDouble4)
    _SET_ARG(CLUByte4)

    void setArg(string function, int arg, CLBuffer * buffer);
    void setArg(string function, int arg, CLBufferGL * buffer);
    void setArg(string function, int arg, CLImageGL * buffer);
    void acquireImageGL(CLImageGL * buffer);
    void releaseImageGL(CLImageGL * buffer);

    ~CLProgram() {
        for (map<string, cl::Kernel*>::iterator ii=functions.begin(); ii!=functions.end(); ii++) {
            delete ii->second;
        }
        functions.clear();
    }
};

class CLBufferGL {
public:
    cl::BufferGL * buffer;
    CLProgram * program;
    GLuint glBuffer;

    CLBufferGL() {
        buffer = NULL;
    }

    CLBufferGL(CLProgram * _program, GLuint bufferGL, MemoryType memType = MEMORY_READ_WRITE) {
        glBuffer = bufferGL;
        buffer = new cl::BufferGL(_program->context->context, static_cast<cl_mem_flags>(memType), glBuffer);
        program = _program;
    }
    ~CLBufferGL() {
        if (buffer != NULL) {
            delete buffer;
            buffer = NULL;
        }
    }
};

class CLImageGL {
public:
    cl::Image2DGL * buffer;
    CLProgram * program;
    GLuint glTex, width, height;
    GLubyte * colorBuffer;

    CLImageGL() {
        buffer = NULL;
        colorBuffer = NULL;
    }

    CLImageGL(CLProgram * _program, GLuint _width, GLuint _height, MemoryType memType = MEMORY_READ_WRITE) {
        width = _width;
        height = _height;
        colorBuffer = new GLubyte[width*height*4];
        glDisable(GL_LIGHTING);
        glEnable(GL_TEXTURE_2D);
        glGenTextures(1, &glTex);
        glBindTexture(GL_TEXTURE_2D, glTex);
        for (size_t i=0; i<width*height; i++) {
            colorBuffer[i*4] = 0;
            colorBuffer[i*4+1] = 0;
            colorBuffer[i*4+2] = 0;
            colorBuffer[i*4+3] = 255;
        }
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, colorBuffer);
        CLInt err;
        buffer = new cl::Image2DGL(_program->context->context, static_cast<cl_mem_flags>(memType), GL_TEXTURE_2D, 0, glTex, &err);
        _program->context->ReportError(err, "CLImageGL(): ");
        program = _program;
    }
    ~CLImageGL() {
        if (buffer != NULL) {
            delete buffer;
            buffer = NULL;
        }
        if (colorBuffer != NULL) {
            delete colorBuffer;
            colorBuffer = NULL;
        }
    }
};

class CLBuffer {
public:
    cl::Buffer * buffer;
    void * data;
    size_t dataSize;
    size_t elSize;
    size_t length;
    CLProgram * program;

    CLBuffer() {
        data = NULL;
        dataSize = 0;
        buffer = NULL;
    }

    CLBuffer(CLProgram * _program, size_t sizeBytes, MemoryType memType = MEMORY_READ_WRITE) {
        elSize = 1;
        length = sizeBytes;
        data = (void*)malloc(dataSize = sizeBytes);
        memset(data, 0, dataSize);
        buffer = new cl::Buffer(_program->context->context, static_cast<cl_mem_flags>(memType), dataSize);
        program = _program;
    }

    CLBuffer(CLProgram * _program, size_t numberElements, size_t elementSize, MemoryType memType = MEMORY_READ_WRITE) {
        elSize = elementSize;
        length = numberElements;
        data = (void*)malloc(dataSize = (numberElements*elementSize));
        memset(data, 0, dataSize);
        buffer = new cl::Buffer(_program->context->context, static_cast<cl_mem_flags>(memType), dataSize);
        program = _program;
    }

    ~CLBuffer() {
        if (data != NULL) {
            free(data);
            data = NULL;
        }
        if (buffer != NULL) {
            delete buffer;
            buffer = NULL;
        }
        dataSize = 0;
    }

    CLFloat * dataFloat () { return (CLFloat*)data; }
    CLDouble * dataDouble () { return (CLDouble*)data; }
    CLUByte * dataUByte () { return (CLUByte*)data; }
    CLFloat2 * dataFloat2 () { return (CLFloat2*)data; }
    CLDouble2 * dataDouble2 () { return (CLDouble2*)data; }
    CLFloat3 * dataFloat3 () { return (CLFloat3*)data; }
    CLDouble3 * dataDouble3 () { return (CLDouble3*)data; }
    CLFloat4 * dataFloat4 () { return (CLFloat4*)data; }
    CLDouble4 * dataDouble4 () { return (CLDouble4*)data; }
    CLUByte4 * dataUByte4 () { return (CLUByte4*)data; }
    CLInt2 * dataInt2 () { return (CLInt2*)data; }
    CLInt3 * dataInt3 () { return (CLInt3*)data; }
    CLInt4 * dataInt4 () { return (CLInt4*)data; }

    void readSync() {
        cl::Event event;
        cl_int err = program->queue.enqueueReadBuffer(*buffer, true, 0, dataSize, data, NULL, &event);
        program->context->ReportError(err, "readSync: ");
        event.wait();
    }

    void writeSync() {
        cl::Event event;
        cl_int err = program->queue.enqueueWriteBuffer(*buffer, true, 0, dataSize, data, NULL, &event);
        program->context->ReportError(err, "writeSync: ");
        event.wait();
    }
};

void CLProgram::setArg(string function, int arg, CLBuffer * buffer) {
    getFunction(function)->setArg(arg, sizeof(cl::Buffer), buffer->buffer);
}

void CLProgram::setArg(string function, int arg, CLBufferGL * buffer) {
    getFunction(function)->setArg(arg, sizeof(cl::BufferGL), buffer->buffer);
}

void CLProgram::setArg(string function, int arg, CLImageGL * buffer) {
    getFunction(function)->setArg(arg, sizeof(cl::Image2DGL), buffer->buffer);
}

void CLProgram::acquireImageGL(CLImageGL * buffer) {
    vector<cl::Memory> mem;
    mem.push_back(*(buffer->buffer));

    cl::Event event;
    CLInt er1 = queue.enqueueAcquireGLObjects(&mem, NULL, &event);
}

void CLProgram::releaseImageGL(CLImageGL * buffer) {
    vector<cl::Memory> mem;
    mem.push_back(*(buffer->buffer));

    cl::Event event;
    
    CLInt er1 = queue.enqueueReleaseGLObjects(&mem, NULL, &event);
    event.wait();
}