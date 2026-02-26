#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <array>
#include <Cell.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "maths/Vec3.hpp"
#include "maths/Mat4.hpp"
#include "maths/Hash.hpp"

#include "Environment.hpp"
#include "Genetic.hpp"
#include "Simulation.hpp"

// ======================= OpenGL helpers ==========================
static GLuint compileShader(const GLenum type, const char* src)
{
    const GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok = 0;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok)
    {
        char log[2048];
        glGetShaderInfoLog(s, sizeof(log), nullptr, log);
        std::cerr << "Shader compile error:\n" << log << "\n";
    }
    return s;
}

static GLuint makeProgram(const char* vs, const char* fs)
{
    const GLuint v = compileShader(GL_VERTEX_SHADER, vs);
    const GLuint f = compileShader(GL_FRAGMENT_SHADER, fs);
    const GLuint p = glCreateProgram();
    glAttachShader(p, v);
    glAttachShader(p, f);
    glLinkProgram(p);
    GLint ok = 0;
    glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok)
    {
        char log[2048];
        glGetProgramInfoLog(p, sizeof(log), nullptr, log);
        std::cerr << "Program link error:\n" << log << "\n";
    }
    glDeleteShader(v);
    glDeleteShader(f);
    return p;
}

// ======================= Render Buffers ==========================
struct RenderPoint
{
    float x, y, z;
    float r, g, b;
    float size;
};

struct RenderLine
{
    float x1, y1, z1;
    float x2, y2, z2;
    float r,  g,  b;
};

// pack lines into same vertex layout as points shader (size unused)
struct LineVertex
{
    float x, y, z;
    float r, g, b;
    float size;
};

int main()
{
    if (!glfwInit())
    {
        std::cerr << "Failed to init GLFW\n";
        return 1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(1400, 900, "EvoLife3D V2 - Dynamic Genome + Multicellular", nullptr, nullptr);
    if (!window)
    {
        std::cerr << "Failed to create window\n";
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress)))
    {
        std::cerr << "Failed to init GLAD\n";
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    auto VS = R"(
        #version 330 core
        layout(location=0) in vec3 aPos;
        layout(location=1) in vec3 aColor;
        layout(location=2) in float aSize;

        uniform mat4 uVP;
        out vec3 vColor;

        void main(){
            gl_Position = uVP * vec4(aPos, 1.0);
            gl_PointSize = aSize;
            vColor = aColor;
        }
    )";

    auto FS_POINTS = R"(
        #version 330 core
        in vec3 vColor;
        out vec4 FragColor;

        void main(){
            vec2 p = gl_PointCoord * 2.0 - 1.0;
            float r2 = dot(p,p);
            if(r2 > 1.0) discard;
            float edge = smoothstep(1.0, 0.5, r2);
            FragColor = vec4(vColor * edge, 1.0);
        }
    )";

    auto FS_LINES = R"(
        #version 330 core
        in vec3 vColor;
        out vec4 FragColor;
        void main(){
            FragColor = vec4(vColor, 0.85);
        }
    )";

    const GLuint progPoints = makeProgram(VS, FS_POINTS);
    const GLuint progLines  = makeProgram(VS, FS_LINES);

    const GLint uVPPoints = glGetUniformLocation(progPoints, "uVP");
    const GLint uVPLines  = glGetUniformLocation(progLines, "uVP");

    // Points VAO/VBO
    GLuint vaoPts = 0, vboPts = 0;
    glGenVertexArrays(1, &vaoPts);
    glGenBuffers(1, &vboPts);

    glBindVertexArray(vaoPts);
    glBindBuffer(GL_ARRAY_BUFFER, vboPts);
    glBufferData(GL_ARRAY_BUFFER, sizeof(RenderPoint) * 1024, nullptr, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(RenderPoint), static_cast<void*>(nullptr));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(RenderPoint), reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(RenderPoint), reinterpret_cast<void*>(6 * sizeof(float)));

    // Lines VAO/VBO (same layout)
    GLuint vaoLines = 0, vboLines = 0;
    glGenVertexArrays(1, &vaoLines);
    glGenBuffers(1, &vboLines);

    glBindVertexArray(vaoLines);
    glBindBuffer(GL_ARRAY_BUFFER, vboLines);
    glBufferData(GL_ARRAY_BUFFER, sizeof(LineVertex) * 1024, nullptr, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(LineVertex), static_cast<void*>(nullptr));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(LineVertex), reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(LineVertex), reinterpret_cast<void*>(6 * sizeof(float)));

    glBindVertexArray(0);

    Sim sim;
    sim.seedInitial(220);

    double       lastTime = glfwGetTime();
    double       accum    = 0.0;
    const double fixedDt  = sim.dt;

    bool  paused       = false;
    bool  showEnv      = true;
    bool  showCurrents = true;
    float camYaw       = 0.8f, camPitch = 0.5f, camDist = 180.f;

    bool pPrev = false, vPrev = false, rPrev = false, nPrev = false, cPrev = false;

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(window, 1);

        // toggles (edge-triggered)
        const bool pNow = (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS);
        if (pNow && !pPrev) paused = !paused;
        pPrev = pNow;

        const bool vNow = (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS);
        if (vNow && !vPrev) showEnv = !showEnv;
        vPrev = vNow;

        const bool cNow = (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS);
        if (cNow && !cPrev) showCurrents = !showCurrents;
        cPrev = cNow;

        const bool rNow = (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS);
        if (rNow && !rPrev)
        {
            sim = Sim{};
            sim.seedInitial(220);
        }
        rPrev = rNow;

        const bool nNow       = (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS);
        const bool singleStep = (nNow && !nPrev && paused);
        nPrev                 = nNow;

        // camera controls
        if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) camYaw -= 0.02f;
        if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) camYaw += 0.02f;
        if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) camPitch = clampf(camPitch + 0.02f, -1.4f, 1.4f);
        if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) camPitch = clampf(camPitch - 0.02f, -1.4f, 1.4f);
        if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) camDist = std::max(40.f, camDist - 2.f);
        if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) camDist = std::min(400.f, camDist + 2.f);

        const double now   = glfwGetTime();
        const double frame = now - lastTime;
        lastTime           = now;
        accum += frame;

        if (!paused)
        {
            int maxSubSteps = 6;
            int steps       = 0;
            while (accum >= fixedDt && steps < maxSubSteps)
            {
                sim.step();
                accum -= fixedDt;
                ++steps;
            }
            if (steps == maxSubSteps && accum > 1.0) accum = 0.0;
        }
        else if (singleStep)
        {
            sim.step();
            accum = 0.0;
        }

        // ---------- Build render lists ----------
        std::vector<RenderPoint> pts;
        std::vector<RenderLine>  lines;
        std::vector<LineVertex>  lineVerts;

        size_t estimatedPts   = 0;
        size_t estimatedLines = 0;
        for (const auto& o : sim.organisms)
        {
            estimatedPts += o.cells.size() + 1;
            estimatedLines += o.springs.size();
        }
        pts.reserve(estimatedPts + 512);
        lines.reserve(estimatedLines);

        // Organisms
        for (const auto& o : sim.organisms)
        {
            const Vec3  baseColor = o.pheno.pigment;
            const float ageT      = clampf(o.age / std::max(1.f, o.pheno.maxAge), 0.f, 1.f);

            // springs
            for (const auto& s : o.springs)
            {
                if (s.a < 0 || s.b < 0 || s.a >= static_cast<int>(o.cells.size()) || s.b >= static_cast<int>(o.cells.
                                                                                                               size()))
                    continue;
                const Vec3 a = o.pos + o.cells[s.a].localPos;
                const Vec3 b = o.pos + o.cells[s.b].localPos;
                const Vec3 col{
                    clampf(0.15f + 0.65f * baseColor.x, 0.f, 1.f),
                    clampf(0.15f + 0.65f * baseColor.y, 0.f, 1.f),
                    clampf(0.15f + 0.65f * baseColor.z, 0.f, 1.f)
                };
                lines.push_back({a.x, a.y, a.z, b.x, b.y, b.z, col.x, col.y, col.z});
            }

            // cells
            for (const auto& c : o.cells)
            {
                const Vec3 wp = o.pos + c.localPos;

                const float n = sim.env.nutrient(wp);
                const float l = sim.env.light(wp);
                const float x = sim.env.toxin(wp);
                const float b = sim.env.biomass(wp);

                const Vec3 col = {
                    clampf(baseColor.x + 0.18f * n + 0.18f * o.pheno.aggression - 0.15f * x, 0.f, 1.f),
                    clampf(baseColor.y + 0.22f * n + 0.10f * b - 0.05f * ageT, 0.f, 1.f),
                    clampf(baseColor.z + 0.24f * l - 0.10f * x, 0.f, 1.f)
                };

                const float size = 3.4f + 2.2f * o.pheno.cellRadius;
                pts.push_back({wp.x, wp.y, wp.z, col.x, col.y, col.z, size});
            }

            // COM marker
            const Vec3 com = Sim::organismCenterOfMass(o);
            pts.push_back({com.x, com.y, com.z, 1.f, 1.f, 1.f, 2.3f});
        }

        // Environment layer: terrain + ocean + currents + field probes
        if (showEnv)
        {
            constexpr int gridN = 26;
            constexpr float gridStep = 3.0f;

            for (int ix = -gridN; ix <= gridN; ++ix)
            {
                for (int iz = -gridN; iz <= gridN; ++iz)
                {
                    const float x = static_cast<float>(ix) * gridStep;
                    const float z = static_cast<float>(iz) * gridStep;
                    const Vec3 sample{x, 0.f, z};

                    const float groundY = sim.env.groundHeight(sample);
                    const float waterY  = sim.env.waterSurfaceHeight(sample);
                    const float shore   = sim.env.shorelineBlend(Vec3{x, waterY, z});
                    const float wetness = clampf((waterY - groundY) / 8.0f, 0.f, 1.f);

                    const float sx = sim.env.groundHeight(Vec3{x + 0.8f, 0.f, z}) - sim.env.groundHeight(
                        Vec3{x - 0.8f, 0.f, z});
                    const float sz = sim.env.groundHeight(Vec3{x, 0.f, z + 0.8f}) - sim.env.groundHeight(
                        Vec3{x, 0.f, z - 0.8f});
                    const float slope = clampf(0.35f * std::sqrt(sx * sx + sz * sz), 0.f, 1.f);

                    // Terrain points
                    const Vec3 landColor{
                        clampf(0.20f + 0.35f * slope + 0.10f * (1.f - wetness), 0.f, 1.f),
                        clampf(0.30f + 0.45f * (1.f - slope) + 0.15f * (1.f - wetness), 0.f, 1.f),
                        clampf(0.16f + 0.18f * (1.f - wetness), 0.f, 1.f)
                    };
                    pts.push_back({x, groundY, z, landColor.x, landColor.y, landColor.z, 1.85f});

                    // Simple terrain contour hints
                    if ((ix % 4 == 0 || iz % 4 == 0) && std::abs(waterY - groundY) < 7.0f)
                    {
                        const float contourY = groundY + 0.08f;
                        const Vec3 contourCol{
                            clampf(0.16f + 0.25f * slope, 0.f, 1.f),
                            clampf(0.20f + 0.20f * (1.f - slope), 0.f, 1.f),
                            0.11f
                        };
                        lines.push_back({x - 0.7f, contourY, z, x + 0.7f, contourY, z, contourCol.x, contourCol.y,
                                         contourCol.z});
                    }

                    // Water surface decoration
                    if (waterY > groundY - 1.3f)
                    {
                        const float wave = 0.5f + 0.5f * std::sin(0.12f * x + 0.09f * z + sim.env.time * 1.2f);
                        const float foam = clampf((waterY - groundY) / 3.5f, 0.f, 1.f) * wave;
                        const Vec3 waterColor{
                            clampf(0.05f + 0.20f * foam, 0.f, 1.f),
                            clampf(0.24f + 0.34f * shore + 0.15f * wave, 0.f, 1.f),
                            clampf(0.40f + 0.46f * shore + 0.10f * wave, 0.f, 1.f)
                        };

                        pts.push_back({x, waterY, z, waterColor.x, waterColor.y, waterColor.z, 2.55f});

                        // shoreline cue
                        if (wetness > 0.10f && wetness < 0.55f)
                        {
                            const Vec3 shoreColor{0.72f, 0.75f, 0.62f};
                            pts.push_back({x, groundY + 0.12f, z, shoreColor.x, shoreColor.y, shoreColor.z, 1.4f});
                        }
                    }
                }
            }

            // Current vectors in water volume
            if (showCurrents)
            {
                for (int ix = -12; ix <= 12; ++ix)
                {
                    for (int iz = -12; iz <= 12; ++iz)
                    {
                        const float x = static_cast<float>(ix) * 5.8f;
                        const float z = static_cast<float>(iz) * 5.8f;

                        for (int iy = 0; iy < 3; ++iy)
                        {
                            const float y = sim.env.fluidLevel - 7.f + static_cast<float>(iy) * 5.0f;
                            const Vec3 p{x, y, z};
                            if (len(p) > sim.worldRadius * 0.95f) continue;

                            const float waterY = sim.env.waterSurfaceHeight(p);
                            const float groundY = sim.env.groundHeight(p);
                            if (p.y > waterY - 0.5f || p.y < groundY + 0.5f) continue;

                            const Vec3 flow = sim.env.fluidVelocity(p);
                            const float flowMag = clampf(len(flow), 0.f, 2.8f);
                            if (flowMag < 0.02f) continue;

                            const Vec3 dir = normalize(flow);
                            const float L = 0.7f + 1.3f * (flowMag / 2.8f);
                            const Vec3 q = p + dir * L;

                            const Vec3 c{
                                clampf(0.08f + 0.18f * flowMag, 0.f, 1.f),
                                clampf(0.45f + 0.18f * flowMag, 0.f, 1.f),
                                clampf(0.62f + 0.24f * flowMag, 0.f, 1.f)
                            };

                            lines.push_back({p.x, p.y, p.z, q.x, q.y, q.z, c.x, c.y, c.z});

                            // arrow head
                            const Vec3 side = normalize(cross(dir, Vec3{0.f, 1.f, 0.f}));
                            const Vec3 h1 = q - dir * 0.35f + side * 0.16f;
                            const Vec3 h2 = q - dir * 0.35f - side * 0.16f;
                            lines.push_back({q.x, q.y, q.z, h1.x, h1.y, h1.z, c.x, c.y, c.z});
                            lines.push_back({q.x, q.y, q.z, h2.x, h2.y, h2.z, c.x, c.y, c.z});
                        }
                    }
                }
            }

            // Sparse volumetric probes (resource fields)
            for (int i = 0; i < 160; i++)
            {
                const Vec3 p{
                    rndf(-sim.worldRadius, sim.worldRadius),
                    rndf(-sim.worldRadius, sim.worldRadius),
                    rndf(-sim.worldRadius, sim.worldRadius)
                };
                if (len(p) > sim.worldRadius) continue;

                const float n = sim.env.nutrient(p);
                const float l = sim.env.light(p);
                const float x = sim.env.toxin(p);
                const float b = sim.env.biomass(p);

                const Vec3 col{
                    clampf(0.10f + 0.50f * x, 0.f, 1.f),
                    clampf(0.10f + 0.52f * n + 0.24f * b, 0.f, 1.f),
                    clampf(0.16f + 0.60f * l, 0.f, 1.f)
                };

                pts.push_back({p.x, p.y, p.z, col.x, col.y, col.z, 1.45f});
            }
        }

        lineVerts.clear();
        lineVerts.reserve(lines.size() * 2);
        for (const auto& [x1, y1, z1, x2, y2, z2, r, g, b] : lines)
        {
            lineVerts.push_back({x1, y1, z1, r, g, b, 1.f});
            lineVerts.push_back({x2, y2, z2, r, g, b, 1.f});
        }

        // Upload
        glBindBuffer(GL_ARRAY_BUFFER, vboPts);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(pts.size() * sizeof(RenderPoint)), pts.data(),
                     GL_DYNAMIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, vboLines);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(lineVerts.size() * sizeof(LineVertex)), lineVerts.data(),
                     GL_DYNAMIC_DRAW);

        // Render
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(0.03f, 0.055f, 0.09f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const Vec3 eye{
            camDist * std::cos(camPitch) * std::cos(camYaw),
            camDist * std::sin(camPitch),
            camDist * std::cos(camPitch) * std::sin(camYaw)
        };

        const Mat4 P = perspective(60.f * 3.1415926f / 180.f,
                                   static_cast<float>(w) / static_cast<float>(std::max(h, 1)),
                                   0.1f, 1000.f);
        const Mat4 V  = lookAt(eye, Vec3{0, 0, 0}, Vec3{0, 1, 0});
        const Mat4 VP = mul(P, V);

        // Lines (springs)
        glUseProgram(progLines);
        glUniformMatrix4fv(uVPLines, 1, GL_FALSE, VP.m);
        glBindVertexArray(vaoLines);
        glLineWidth(1.0f);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(lineVerts.size()));

        // Points (cells + env probes)
        glUseProgram(progPoints);
        glUniformMatrix4fv(uVPPoints, 1, GL_FALSE, VP.m);
        glBindVertexArray(vaoPts);
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(pts.size()));

        glfwSwapBuffers(window);

        // Console stats every ~1 sec
        static double statTimer = 0.0;
        statTimer += frame;
        if (statTimer > 1.0)
        {
            statTimer           = 0.0;
            size_t totalCells   = 0;
            size_t totalSprings = 0;
            for (const auto& o : sim.organisms)
            {
                totalCells += o.cells.size();
                totalSprings += o.springs.size();
            }

            std::cout << "Organisms: " << sim.organisms.size()
                << " | Cells: " << totalCells
                << " | Springs: " << totalSprings
                << " | BiomassPulses: " << sim.env.biomassPulses.size()
                << " | Time: " << sim.env.time
                << (paused ? " [PAUSED]" : "")
                << (showEnv ? " [ENV ON]" : " [ENV OFF]")
                << (showCurrents ? " [CUR ON]" : " [CUR OFF]")
                << "\n";
        }
    }

    glDeleteBuffers(1, &vboPts);
    glDeleteVertexArrays(1, &vaoPts);
    glDeleteBuffers(1, &vboLines);
    glDeleteVertexArrays(1, &vaoLines);
    glDeleteProgram(progPoints);
    glDeleteProgram(progLines);

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
