#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "evo/Math.hpp"
#include "evo/Renderer.hpp"
#include "evo/Sim.hpp"

using namespace evo;

static void glfwErrorCallback(const int code, const char *msg)
{
    std::cerr << "GLFW error (" << code << "): " << (msg ? msg : "") << "\n";
}

int main()
{
    glfwSetErrorCallback(glfwErrorCallback);

    if (!glfwInit())
    {
        std::cerr << "Failed to init GLFW\n";
        return 1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow *window = glfwCreateWindow(1400, 900, "Genetic-Sim V3 (Instanced Spheres)", nullptr, nullptr);
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

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    Renderer renderer;
    if (!renderer.init())
    {
        std::cerr << "Failed to init renderer\n";
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    Sim sim(1234567ULL);
    sim.seedInitial(220);

    bool  paused  = false;
    bool  showEnv = true;
    float camYaw  = 0.8f, camPitch = 0.45f, camDist = 185.f;

    double       lastTime = glfwGetTime();
    double       accum    = 0.0;
    const double fixedDt  = sim.dt;

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, 1);

        // toggles
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
            paused = true;
        if (glfwGetKey(window, GLFW_KEY_ENTER) == GLFW_PRESS)
            paused = false;

        if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS)
            renderer.setShowLines(true);
        if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS)
            renderer.setShowLines(false);

        if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS)
            showEnv = true;
        if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS)
            showEnv = false;

        // camera
        if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
            camYaw -= 0.02f;
        if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
            camYaw += 0.02f;
        if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
            camPitch = clampf(camPitch + 0.02f, -1.35f, 1.35f);
        if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
            camPitch = clampf(camPitch - 0.02f, -1.35f, 1.35f);
        if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
            camDist = std::max(45.f, camDist - 2.f);
        if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
            camDist = std::min(520.f, camDist + 2.f);

        // fixed step
        const double now   = glfwGetTime();
        const double frame = now - lastTime;
        lastTime           = now;
        accum += frame;

        if (!paused)
        {
            int           steps       = 0;
            constexpr int maxSubSteps = 5;
            while (accum >= fixedDt && steps < maxSubSteps)
            {
                sim.step();
                accum -= fixedDt;
                steps++;
            }
        }

        // Build instance list
        std::vector<RenderPoint> inst;
        size_t                   estimated = 0;
        for (const auto &o : sim.organisms())
            estimated += o.cells.size();
        inst.reserve(estimated + 512);

        // optional debug lines (springs)
        std::vector<LineVertex> lines;
        if (renderer.showLines())
        {
            lines.reserve(estimated * 2);
        }

        for (const auto &o : sim.organisms())
        {
            const Vec3  base = o.pheno.pigment;
            const float ageT = clampf(o.age / o.pheno.maxAge, 0.f, 1.f);

            // cells
            for (int i = 0; i < static_cast<int>(o.cells.size()); ++i)
            {
                const auto &c  = o.cells[static_cast<size_t>(i)];
                const Vec3  wp = o.pos + c.localPos;

                const float n = sim.env().nutrient(wp);
                const float l = sim.env().light(wp);
                const float x = sim.env().toxin(wp);
                const float b = sim.env().biomass(wp);

                Vec3 col{
                    clampf(base.x + 0.18f * n + 0.22f * o.pheno.aggression - 0.14f * x + 0.06f * b, 0.f, 1.f),
                    clampf(base.y + 0.25f * n - 0.05f * ageT + 0.05f * b, 0.f, 1.f),
                    clampf(base.z + 0.25f * l - 0.10f * x, 0.f, 1.f),
                };

                const float radius = 0.55f * o.pheno.cellRadius; // world scale
                inst.push_back({wp.x, wp.y, wp.z, col.x, col.y, col.z, radius});

                if (renderer.showLines())
                {
                    // add springs as lines
                    // (we draw them once per spring below, not per cell)
                }
            }

            if (renderer.showLines())
            {
                for (const auto &s : o.springs)
                {
                    if (s.a < 0 || s.b < 0 || s.a >= static_cast<int>(o.cells.size()) || s.b >= static_cast<
                            int>(o.cells.size()))
                        continue;
                    const Vec3 A = o.pos + o.cells[static_cast<size_t>(s.a)].localPos;
                    const Vec3 B = o.pos + o.cells[static_cast<size_t>(s.b)].localPos;

                    constexpr Vec3 lc{0.9f, 0.9f, 0.95f};
                    lines.push_back({A.x, A.y, A.z, lc.x, lc.y, lc.z, 1.f});
                    lines.push_back({B.x, B.y, B.z, lc.x, lc.y, lc.z, 1.f});
                }
            }
        }

        // Environment probes (optional)
        if (showEnv)
        {
            for (int i = 0; i < 260; ++i)
            {
                Vec3 p{
                    static_cast<float>(std::sin(0.17 * i + sim.env().time * 0.11f)) * sim.worldRadius * 0.85f,
                    static_cast<float>(std::sin(0.13 * i + sim.env().time * 0.09f)) * sim.worldRadius * 0.85f,
                    static_cast<float>(std::cos(0.19 * i + sim.env().time * 0.10f)) * sim.worldRadius * 0.85f,
                };

                const float n = sim.env().nutrient(p);
                const float l = sim.env().light(p);
                const float x = sim.env().toxin(p);
                const float b = sim.env().biomass(p);

                Vec3 col{
                    clampf(0.20f + 0.55f * x + 0.15f * b, 0.f, 1.f), // red = toxin + biomass
                    clampf(0.15f + 0.60f * n + 0.10f * b, 0.f, 1.f), // green = nutrient + biomass
                    clampf(0.15f + 0.60f * l, 0.f, 1.f),             // blue = light
                };

                inst.push_back({p.x, p.y, p.z, col.x, col.y, col.z, 0.25f});
            }
        }

        // Upload to GPU
        renderer.uploadInstances(inst);
        if (renderer.showLines())
            renderer.uploadLines(lines);

        // Camera matrices
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        renderer.resize(w, h);

        const Vec3 eye{
            camDist * std::cos(camPitch) * std::cos(camYaw),
            camDist * std::sin(camPitch),
            camDist * std::cos(camPitch) * std::sin(camYaw),
        };

        const Mat4 P = perspective(60.f * 3.1415926535f / 180.f,
                                   static_cast<float>(w) / static_cast<float>(std::max(h, 1)), 0.1f,
                                   1400.f);
        const Mat4 V  = lookAt(eye, Vec3{0, 0, 0}, Vec3{0, 1, 0});
        const Mat4 VP = mul(P, V);
        renderer.setCameraVP(VP);

        // Render
        glClearColor(0.02f, 0.02f, 0.03f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        renderer.drawSpheres();
        renderer.drawLines();

        glfwSwapBuffers(window);

        // Stats
        static double statTimer = 0.0;
        statTimer += frame;
        if (statTimer > 1.0)
        {
            statTimer         = 0.0;
            size_t totalCells = 0;
            for (const auto &o : sim.organisms())
                totalCells += o.cells.size();

            std::cout << "Organisms: " << sim.organisms().size() << " | Cells: " << totalCells
                << " | EnvT: " << sim.env().time << (paused ? " [PAUSED]" : "")
                << " | Lines: " << (renderer.showLines() ? "ON" : "OFF")
                << " | EnvViz: " << (showEnv ? "ON" : "OFF") << "\n";
        }
    }

    renderer.shutdown();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}