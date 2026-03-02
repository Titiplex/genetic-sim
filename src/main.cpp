#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <imgui.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>

#include "ui/Picking.hpp"

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

    // ===== ImGui init =====
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    // ===== Renderer =====
    Renderer renderer;
    if (!renderer.init())
    {
        std::cerr << "Failed to init renderer\n";
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    // ===== Sim =====
    Sim sim(1234567ULL);
    sim.seedInitial(220);

    // ===== UI / runtime flags =====
    bool paused   = false;
    bool stepOnce = false;

    bool showUI  = true;
    bool showEnv = true;
    bool showCom = true;

    float camYaw   = 0.8f;
    float camPitch = 0.45f;
    float camDist  = 185.f;

    int      initialPop = 260;
    int      targetPop  = sim.targetPopulation;
    float    simDt      = sim.dt;
    uint64_t seed       = 1234567ULL;

    int selectedOrg = -1; // index in sim.organisms()

    // rolling history
    static constexpr int HISTORY = 240;
    static float         histPop[HISTORY]{};
    static float         histCells[HISTORY]{};
    static int           histHead = 0;

    auto pushHist = [&](const float pop, const float cells)
    {
        histPop[histHead]   = pop;
        histCells[histHead] = cells;
        histHead            = (histHead + 1) % HISTORY;
    };

    // edge-triggered input helpers
    auto keyPressedOnce = [&](const int key) -> bool
    {
        static bool prev[512]{};
        const bool  cur = glfwGetKey(window, key) == GLFW_PRESS;
        const bool  out = cur && !prev[key];
        prev[key]       = cur;
        return out;
    };

    auto mouseLeftPressedOnce = [&]() -> bool
    {
        static bool prev = false;
        const bool  cur  = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        const bool  out  = cur && !prev;
        prev             = cur;
        return out;
    };

    // fixed step timing
    double lastTime = glfwGetTime();
    double accum    = 0.0;
    double fixedDt  = sim.dt;

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        // ---- global exit
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, 1);

        // ---- edge-triggered toggles
        if (keyPressedOnce(GLFW_KEY_F1))
            showUI = !showUI;
        if (keyPressedOnce(GLFW_KEY_SPACE))
            paused = !paused;
        if (keyPressedOnce(GLFW_KEY_S))
            stepOnce = true;

        if (keyPressedOnce(GLFW_KEY_L))
            renderer.setShowLines(true);
        if (keyPressedOnce(GLFW_KEY_K))
            renderer.setShowLines(false);

        if (keyPressedOnce(GLFW_KEY_O))
            showEnv = true;
        if (keyPressedOnce(GLFW_KEY_P))
            showEnv = false;

        // ---- camera continuous
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

        // ---- fixed step
        const double now   = glfwGetTime();
        const double frame = now - lastTime;
        lastTime           = now;
        accum += frame;

        if (!paused || stepOnce)
        {
            constexpr int maxSubSteps = 4;
            int           steps       = 0;

            if (stepOnce)
            {
                accum += fixedDt;
                stepOnce = false;
            }

            while (accum >= fixedDt && steps < maxSubSteps)
            {
                sim.step();
                accum -= fixedDt;
                steps++;
            }
        }

        // ---- build instances
        std::vector<RenderPoint> inst;
        size_t                   estimated = 0;
        for (const auto &o : sim.organisms())
            estimated += o.cells.size();
        inst.reserve(estimated + 1024);

        std::vector<LineVertex> lines;
        if (renderer.showLines())
            lines.reserve(estimated * 2);

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

                const auto &mode = o.pheno.modes[static_cast<size_t>(c.modeId)];
                col.x            = clampf(col.x * mode.tint.x, 0.f, 1.f);
                col.y            = clampf(col.y * mode.tint.y, 0.f, 1.f);
                col.z            = clampf(col.z * mode.tint.z, 0.f, 1.f);

                const float radius = 0.55f * o.pheno.cellRadius;
                inst.push_back({wp.x, wp.y, wp.z, col.x, col.y, col.z, radius});
            }

            // COM marker
            if (showCom)
            {
                const Vec3 com = Sim::centerOfMass(o);
                inst.push_back({com.x, com.y, com.z, 1.f, 1.f, 1.f, 0.18f});
            }

            // springs
            if (renderer.showLines())
            {
                constexpr Vec3 lc{0.9f, 0.9f, 0.95f};
                for (const auto &s : o.springs)
                {
                    if (s.a < 0 || s.b < 0 || s.a >= static_cast<int>(o.cells.size()) || s.b >= static_cast<
                            int>(o.cells.size()))
                        continue;

                    const Vec3 A = o.pos + o.cells[static_cast<size_t>(s.a)].localPos;
                    const Vec3 B = o.pos + o.cells[static_cast<size_t>(s.b)].localPos;
                    lines.push_back({A.x, A.y, A.z, lc.x, lc.y, lc.z, 1.f});
                    lines.push_back({B.x, B.y, B.z, lc.x, lc.y, lc.z, 1.f});
                }
            }
        }

        // environment probes
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
                    clampf(0.20f + 0.55f * x + 0.15f * b, 0.f, 1.f),
                    clampf(0.15f + 0.60f * n + 0.10f * b, 0.f, 1.f),
                    clampf(0.15f + 0.60f * l, 0.f, 1.f),
                };

                inst.push_back({p.x, p.y, p.z, col.x, col.y, col.z, 0.25f});
            }
        }

        // upload
        renderer.uploadInstances(inst);
        if (renderer.showLines())
            renderer.uploadLines(lines);

        // camera matrices
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        renderer.resize(w, h);

        const Vec3 eye{
            camDist * std::cos(camPitch) * std::cos(camYaw),
            camDist * std::sin(camPitch),
            camDist * std::cos(camPitch) * std::sin(camYaw),
        };

        const Mat4 P = perspective(60.f * 3.1415926535f / 180.f,
                                   static_cast<float>(w) / static_cast<float>(std::max(h, 1)), 0.1f, 1400.f);
        const Mat4 V  = lookAt(eye, Vec3{0, 0, 0}, Vec3{0, 1, 0});
        const Mat4 VP = mul(P, V);
        renderer.setCameraVP(VP);

        // picking (click once)
        if (mouseLeftPressedOnce() && !ImGui::GetIO().WantCaptureMouse)
        {
            double mx, my;
            glfwGetCursorPos(window, &mx, &my);

            const auto &orgs    = sim.organisms();
            float       best    = 1e9f;
            int         bestIdx = -1;

            for (int i = 0; i < static_cast<int>(orgs.size()); ++i)
            {
                const auto &o   = orgs[static_cast<size_t>(i)];
                const Vec3  com = Sim::centerOfMass(o);

                Vec2 sc;
                if (float depth; !worldToScreen(com, VP, w, h, sc, depth))
                    continue;

                const float dx = static_cast<float>(mx) - sc.x;
                const float dy = static_cast<float>(my) - sc.y;

                if (const float d2 = dx * dx + dy * dy; d2 < best && d2 < 18.f * 18.f)
                {
                    best    = d2;
                    bestIdx = i;
                }
            }
            selectedOrg = bestIdx;
        }

        // ---- Render scene
        glClearColor(0.02f, 0.02f, 0.03f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        renderer.drawSpheres();
        renderer.drawLines();

        // ---- UI
        if (showUI)
        {
            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();

            ImGui::Begin("EvoLife3D");

            ImGui::Text("F1 toggle UI | SPACE pause | S step | click select");
            ImGui::Separator();

            if (ImGui::Button(paused ? "Run" : "Pause"))
                paused = !paused;
            ImGui::SameLine();
            if (ImGui::Button("Step"))
                stepOnce = true;
            ImGui::SameLine();
            if (ImGui::Button("Reset"))
            {
                sim.reset(seed, initialPop);
                selectedOrg = -1;
                for (int i = 0; i < HISTORY; ++i)
                {
                    histPop[i]   = 0.f;
                    histCells[i] = 0.f;
                }
                histHead = 0;
            }

            ImGui::InputScalar("Seed", ImGuiDataType_U64, &seed);
            ImGui::SliderInt("Initial Pop", &initialPop, 50, 1200);

            ImGui::SliderInt("Target Pop", &targetPop, 50, 4000);
            sim.targetPopulation = targetPop;

            ImGui::SliderFloat("dt", &simDt, 0.005f, 0.08f, "%.3f");
            sim.dt  = simDt;
            fixedDt = sim.dt;

            ImGui::Checkbox("Show env points", &showEnv);
            ImGui::Checkbox("Show COM markers", &showCom);

            bool showLines = renderer.showLines();
            if (ImGui::Checkbox("Show springs (lines)", &showLines))
                renderer.setShowLines(showLines);

            ImGui::Separator();

            const auto [organisms, totalCells, avgCells, modeCounts] = sim.computeStats();
            pushHist(static_cast<float>(organisms), static_cast<float>(totalCells));

            ImGui::Text("Organisms: %d", organisms);
            ImGui::Text("Total cells: %d (avg %.2f)", totalCells, avgCells);

            static float popPlot[HISTORY];
            static float cellPlot[HISTORY];
            for (int i = 0; i < HISTORY; ++i)
            {
                const int idx = (histHead + i) % HISTORY;
                popPlot[i]    = histPop[idx];
                cellPlot[i]   = histCells[idx];
            }

            ImGui::PlotLines("Pop", popPlot, HISTORY, 0, nullptr, 0.f, 4000.f, ImVec2(0, 70));
            ImGui::PlotLines("Cells", cellPlot, HISTORY, 0, nullptr, 0.f, 200000.f, ImVec2(0, 70));

            ImGui::Separator();
            ImGui::Text("Mode distribution:");
            for (int k = 0; k < Phenotype::K; ++k)
                ImGui::Text("  mode %d: %d", k, modeCounts[static_cast<size_t>(k)]);

            ImGui::Separator();
            if (selectedOrg >= 0 && selectedOrg < static_cast<int>(sim.organisms().size()))
            {
                const auto &o = sim.organisms()[static_cast<size_t>(selectedOrg)];
                ImGui::Text("Selected: #%llu", static_cast<unsigned long long>(o.id));
                ImGui::Text("Cells: %d  Energy: %.2f  Age: %.2f/%.2f",
                            static_cast<int>(o.cells.size()), o.energy, o.age, o.pheno.maxAge);

                ImGui::Text("Phenotype:");
                ImGui::Text("  nutrientAffinity: %.2f", o.pheno.nutrientAffinity);
                ImGui::Text("  photoAffinity: %.2f", o.pheno.photoAffinity);
                ImGui::Text("  toxinResistance: %.2f", o.pheno.toxinResistance);
                ImGui::Text("  aggression: %.2f", o.pheno.aggression);
                ImGui::Text("  growthBias: %.2f", o.pheno.growthBias);
                ImGui::Text("  gateSharpness: %.2f", o.pheno.gateSharpness);

                if (ImGui::TreeNode("Genome"))
                {
                    for (size_t i = 0; i < o.genome.genes.size(); ++i)
                        ImGui::Text("gene[%zu] = %s", i, o.genome.genes[i].allele.c_str());
                    ImGui::TreePop();
                }

                if (ImGui::TreeNode("Modes"))
                {
                    for (int k = 0; k < Phenotype::K; ++k)
                    {
                        const auto &m = o.pheno.modes[static_cast<size_t>(k)];
                        ImGui::Text(
                            "Mode %d: UptN %.2f  Photo %.2f  DigB %.2f  ResX %.2f  Atk %.2f  Adh %.2f  Maint %.2f",
                            k, m.gUptakeN, m.gPhoto, m.gDigestB, m.gResistX, m.gAttack, m.gAdhesion,
                            m.gMaint);
                    }
                    ImGui::TreePop();
                }
            }
            else
            {
                ImGui::Text("Selected: none (click an organism)");
            }

            ImGui::End();

            ImGui::Render();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        }

        // swap LAST
        glfwSwapBuffers(window);

        // console stats ~1sec
        static double statTimer = 0.0;
        statTimer += frame;
        if (statTimer > 1.0)
        {
            statTimer         = 0.0;
            size_t totalCells = 0;
            for (const auto &o : sim.organisms())
                totalCells += o.cells.size();

            std::cout
                << "Organisms: " << sim.organisms().size()
                << " | Cells: " << totalCells
                << " | EnvT: " << sim.env().time
                << (paused ? " [PAUSED]" : "")
                << " | Lines: " << (renderer.showLines() ? "ON" : "OFF")
                << " | EnvViz: " << (showEnv ? "ON" : "OFF")
                << "\n";
        }
    }

    renderer.shutdown();

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}