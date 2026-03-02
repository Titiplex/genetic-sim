#pragma once
#include "evo/Math.hpp" // là où tu as Vec3/Mat4/mul etc

struct Vec2
{
    float x = 0, y = 0;
};

static inline bool worldToScreen(const evo::Vec3 &wpos, const evo::Mat4 &VP, const int w, const int h,
                                 Vec2 &           out, float &           outDepth)
{
    // OpenGL clip coords: clip = VP * vec4(wpos,1)
    const float  x = wpos.x, y = wpos.y, z = wpos.z;
    const float *m = VP.m;

    // Column-major: result = M * v
    const float cx = m[0] * x + m[4] * y + m[8] * z + m[12] * 1.f;
    const float cy = m[1] * x + m[5] * y + m[9] * z + m[13] * 1.f;
    const float cz = m[2] * x + m[6] * y + m[10] * z + m[14] * 1.f;
    const float cw = m[3] * x + m[7] * y + m[11] * z + m[15] * 1.f;

    if (cw <= 1e-6f)
        return false;

    const float ndcX = cx / cw;
    const float ndcY = cy / cw;
    const float ndcZ = cz / cw;

    // behind / outside clip can still be shown; we threshold lightly
    out.x    = (ndcX * 0.5f + 0.5f) * static_cast<float>(w);
    out.y    = (1.f - (ndcY * 0.5f + 0.5f)) * static_cast<float>(h); // flip for screen coords
    outDepth = ndcZ;
    return true;
}