#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <array>
#include <Cell.hpp>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <deque>
#include <iomanip>
#include <random>
#include <sstream>
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

enum class AppScreen : uint8_t
{
    MainMenu,
    Simulation,
    Options
};

static Vec3 biomeColor(const Biome biome)
{
    switch (biome)
    {
        case Biome::OpenOcean: return {0.05f, 0.28f, 0.62f};
        case Biome::ReefShelf: return {0.10f, 0.44f, 0.68f};
        case Biome::KelpLagoon: return {0.12f, 0.52f, 0.44f};
        case Biome::Estuary: return {0.24f, 0.48f, 0.36f};
        case Biome::Mangrove: return {0.22f, 0.40f, 0.22f};
        case Biome::TidalFlat: return {0.44f, 0.42f, 0.28f};
        case Biome::CoastalForest: return {0.17f, 0.50f, 0.19f};
        case Biome::InlandWetland: return {0.18f, 0.56f, 0.34f};
        case Biome::Savanna: return {0.58f, 0.52f, 0.26f};
        case Biome::Alpine: return {0.62f, 0.64f, 0.66f};
    }
    return {0.4f, 0.4f, 0.4f};
}

struct AppSettings
{
    int   initialPopulation = 220;
    float startTimeScale    = 1.0f;
    bool  startPaused       = false;
    bool  showEnv           = true;
    bool  showCurrents      = true;
    bool  cinematicCam      = true;
    bool  showTrails        = true;
};

static void resetSimulation(Sim& sim, const AppSettings& settings, float& timeScale, bool& paused)
{
    sim = Sim{};
    sim.seedInitial(settings.initialPopulation);
    timeScale = settings.startTimeScale;
    paused    = settings.startPaused;
}

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
        out float vViewZ;

        void main(){
            gl_Position = uVP * vec4(aPos, 1.0);
            gl_PointSize = aSize;
            vColor = aColor;
            vViewZ = gl_Position.z / max(0.0001, gl_Position.w);
        }
    )";

    auto FS_POINTS = R"(
        #version 330 core
        in vec3 vColor;
        in float vViewZ;
        out vec4 FragColor;

        uniform float uTime;
        uniform float uFogDensity;
        uniform float uGlow;

        void main(){
            vec2 p = gl_PointCoord * 2.0 - 1.0;
            float r2 = dot(p,p);
            if(r2 > 1.0) discard;
            float radial = 1.0 - clamp(r2, 0.0, 1.0);
            float core = smoothstep(0.0, 1.0, radial);
            float halo = smoothstep(0.0, 0.65, radial);
            float pulse = 0.92 + 0.08 * sin(uTime * 1.7 + (vColor.r + vColor.g + vColor.b) * 5.0);
            float fog = exp(-uFogDensity * (vViewZ * vViewZ + 0.25));
            vec3 lit = vColor * (0.35 + 0.85 * core) + vec3(0.18, 0.22, 0.26) * halo * uGlow;
            FragColor = vec4(lit * fog * pulse, clamp(0.4 + halo * 0.65, 0.0, 1.0));
        }
    )";

    auto FS_LINES = R"(
        #version 330 core
        in vec3 vColor;
        in float vViewZ;
        out vec4 FragColor;
        uniform float uFogDensity;
        void main(){
            float fog = exp(-uFogDensity * (vViewZ * vViewZ + 0.15));
            FragColor = vec4(vColor * fog, 0.72);
        }
    )";

    const GLuint progPoints = makeProgram(VS, FS_POINTS);
    const GLuint progLines  = makeProgram(VS, FS_LINES);

    const GLint uVPPoints = glGetUniformLocation(progPoints, "uVP");
    const GLint uVPLines  = glGetUniformLocation(progLines, "uVP");
    const GLint uTimePoints = glGetUniformLocation(progPoints, "uTime");
    const GLint uFogPoints = glGetUniformLocation(progPoints, "uFogDensity");
    const GLint uFogLines = glGetUniformLocation(progLines, "uFogDensity");
    const GLint uGlowPoints = glGetUniformLocation(progPoints, "uGlow");

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

    Sim         sim;
    AppSettings settings;
    float       timeScale = settings.startTimeScale;
    bool        paused    = settings.startPaused;
    resetSimulation(sim, settings, timeScale, paused);

    double       lastTime = glfwGetTime();
    double       accum    = 0.0;
    const double fixedDt  = sim.dt;

    bool  showEnv      = settings.showEnv;
    bool  showCurrents = settings.showCurrents;
    bool  cinematicCam = settings.cinematicCam;
    bool  showTrails   = settings.showTrails;
    AppScreen appScreen = AppScreen::MainMenu;
    AppScreen optionsReturnScreen = AppScreen::MainMenu;
    bool pausedBeforeOptions = paused;
    float camYaw       = 0.8f, camPitch = 0.5f, camDist = 180.f;
    float fogDensity   = 0.65f;

    std::unordered_map<uint64_t, std::deque<Vec3>> trails;

    bool pPrev = false, vPrev = false, rPrev = false, nPrev = false, cPrev = false, fPrev = false, tPrev = false;
    bool enterPrev = false, mPrev = false, oPrev = false;
    bool iPrev = false, kPrev = false, gPrev = false;
    bool onePrev = false, twoPrev = false;

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) glfwSetWindowShouldClose(window, 1);

        // Global app navigation
        const bool enterNow = (glfwGetKey(window, GLFW_KEY_ENTER) == GLFW_PRESS);
        if (enterNow && !enterPrev && appScreen == AppScreen::MainMenu)
        {
            appScreen = AppScreen::Simulation;
            resetSimulation(sim, settings, timeScale, paused);
            showEnv      = settings.showEnv;
            showCurrents = settings.showCurrents;
            cinematicCam = settings.cinematicCam;
            showTrails   = settings.showTrails;
            trails.clear();
            accum = 0.0;
        }
        enterPrev = enterNow;

        const bool oNow = (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS);
        if (oNow && !oPrev)
        {
            if (appScreen == AppScreen::MainMenu || appScreen == AppScreen::Simulation)
            {
                optionsReturnScreen = appScreen;
                pausedBeforeOptions = paused;
                appScreen = AppScreen::Options;
                paused = true;
            }
            else if (appScreen == AppScreen::Options)
            {
                appScreen = optionsReturnScreen;
                if (appScreen == AppScreen::Simulation) paused = pausedBeforeOptions;
            }
        }
        oPrev = oNow;

        const bool mNow = (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS);
        if (mNow && !mPrev)
        {
            if (appScreen != AppScreen::MainMenu)
            {
                appScreen = AppScreen::MainMenu;
                paused = true;
            }
        }
        mPrev = mNow;

        // options menu controls
        if (appScreen == AppScreen::Options)
        {
            const bool iNow = (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS);
            if (iNow && !iPrev) settings.initialPopulation = std::min(500, settings.initialPopulation + 10);
            iPrev = iNow;

            const bool kNow = (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS);
            if (kNow && !kPrev) settings.initialPopulation = std::max(20, settings.initialPopulation - 10);
            kPrev = kNow;

            const bool gNow = (glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS);
            if (gNow && !gPrev)
            {
                settings.showEnv = !settings.showEnv;
                showEnv = settings.showEnv;
            }
            gPrev = gNow;

            const bool cNowOpt = (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS);
            if (cNowOpt && !cPrev)
            {
                settings.showCurrents = !settings.showCurrents;
                showCurrents = settings.showCurrents;
            }
            cPrev = cNowOpt;

            const bool fNowOpt = (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS);
            if (fNowOpt && !fPrev)
            {
                settings.cinematicCam = !settings.cinematicCam;
                cinematicCam = settings.cinematicCam;
            }
            fPrev = fNowOpt;

            const bool tNowOpt = (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS);
            if (tNowOpt && !tPrev)
            {
                settings.showTrails = !settings.showTrails;
                showTrails = settings.showTrails;
            }
            tPrev = tNowOpt;

            const bool pNowOpt = (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS);
            if (pNowOpt && !pPrev)
            {
                settings.startPaused = !settings.startPaused;
            }
            pPrev = pNowOpt;

            const bool oneNow = (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS);
            if (oneNow && !onePrev) settings.startTimeScale = std::max(0.2f, settings.startTimeScale - 0.2f);
            onePrev = oneNow;

            const bool twoNow = (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS);
            if (twoNow && !twoPrev) settings.startTimeScale = std::min(4.0f, settings.startTimeScale + 0.2f);
            twoPrev = twoNow;
        }
        else
        {
            iPrev = kPrev = gPrev = onePrev = twoPrev = false;
        }

        const bool simulationRunning = (appScreen == AppScreen::Simulation);

        // toggles (edge-triggered)
        const bool pNow = (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS);
        if (simulationRunning && pNow && !pPrev) paused = !paused;
        pPrev = pNow;

        const bool vNow = (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS);
        if (simulationRunning && vNow && !vPrev) showEnv = !showEnv;
        vPrev = vNow;

        const bool cNow = (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS);
        if (simulationRunning && cNow && !cPrev) showCurrents = !showCurrents;
        cPrev = cNow;

        const bool rNow = (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS);
        if (simulationRunning && rNow && !rPrev)
        {
            resetSimulation(sim, settings, timeScale, paused);
            showEnv      = settings.showEnv;
            showCurrents = settings.showCurrents;
            cinematicCam = settings.cinematicCam;
            showTrails   = settings.showTrails;
            trails.clear();
            accum = 0.0;
        }
        rPrev = rNow;

        const bool nNow       = (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS);
        const bool singleStep = (simulationRunning && nNow && !nPrev && paused);
        nPrev                 = nNow;

        const bool fNow = (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS);
        if (simulationRunning && fNow && !fPrev) cinematicCam = !cinematicCam;
        fPrev = fNow;

        const bool tNow = (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS);
        if (simulationRunning && tNow && !tPrev) showTrails = !showTrails;
        tPrev = tNow;

        if (simulationRunning && (glfwGetKey(window, GLFW_KEY_KP_ADD) == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_EQUAL) == GLFW_PRESS)) timeScale = std::min(6.f, timeScale + 0.015f);
        if (simulationRunning && (glfwGetKey(window, GLFW_KEY_KP_SUBTRACT) == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_MINUS) == GLFW_PRESS)) timeScale = std::max(0.1f, timeScale - 0.015f);
        if (simulationRunning && glfwGetKey(window, GLFW_KEY_LEFT_BRACKET) == GLFW_PRESS) fogDensity = std::max(0.08f, fogDensity - 0.008f);
        if (simulationRunning && glfwGetKey(window, GLFW_KEY_RIGHT_BRACKET) == GLFW_PRESS) fogDensity = std::min(1.7f, fogDensity + 0.008f);

        // camera controls
        if (simulationRunning && glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) camYaw -= 0.02f;
        if (simulationRunning && glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) camYaw += 0.02f;
        if (simulationRunning && glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) camPitch = clampf(camPitch + 0.02f, -1.4f, 1.4f);
        if (simulationRunning && glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) camPitch = clampf(camPitch - 0.02f, -1.4f, 1.4f);
        if (simulationRunning && glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) camDist = std::max(40.f, camDist - 2.f);
        if (simulationRunning && glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) camDist = std::min(400.f, camDist + 2.f);

        const double now   = glfwGetTime();
        const double frame = now - lastTime;
        lastTime           = now;
        accum += frame * timeScale;

        if (simulationRunning && !paused)
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
        else if (simulationRunning && singleStep)
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
            estimatedPts += (o.macroMode ? 18 : o.cells.size()) + 2;
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

            // cells / macro representation
            const Vec3 com = Sim::organismCenterOfMass(o);
            if (o.macroMode)
            {
                const float marine = clampf((sim.env.waterSurfaceHeight(com) - com.y) / 8.f, 0.f, 1.f);
                const Vec3 shell = {
                    clampf(baseColor.x + 0.22f * marine + 0.12f * o.terrestrialAffinity, 0.f, 1.f),
                    clampf(baseColor.y + 0.12f * marine, 0.f, 1.f),
                    clampf(baseColor.z + 0.28f * (1.f - o.terrestrialAffinity), 0.f, 1.f)
                };
                pts.push_back({com.x, com.y, com.z, shell.x, shell.y, shell.z, 8.0f + 5.2f * o.macroScale});
                for (int k = 0; k < 12; ++k)
                {
                    const Vec3 ring = normalize(hashToVec3(o.id + static_cast<uint64_t>(k * 17 + 900)));
                    const Vec3 rp = com + ring * (1.4f + 1.1f * o.macroScale);
                    pts.push_back({rp.x, rp.y, rp.z, shell.x * 0.85f, shell.y * 0.9f, shell.z, 3.2f + 0.9f * o.macroScale});
                    lines.push_back({com.x, com.y, com.z, rp.x, rp.y, rp.z, shell.x * 0.55f, shell.y * 0.65f, shell.z * 0.8f});
                }
            }
            else
            {
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
            }

            // COM marker
            pts.push_back({com.x, com.y, com.z, 1.f, 1.f, 1.f, 2.3f});

            auto& trail = trails[o.id];
            trail.push_back(com);
            if (trail.size() > 22) trail.pop_front();
            if (showTrails && trail.size() > 2)
            {
                for (size_t k = 1; k < trail.size(); ++k)
                {
                    const float t = static_cast<float>(k) / static_cast<float>(trail.size() - 1);
                    const Vec3 tc = {
                        clampf(baseColor.x * (0.2f + 0.6f * t), 0.f, 1.f),
                        clampf(baseColor.y * (0.2f + 0.6f * t), 0.f, 1.f),
                        clampf(baseColor.z * (0.2f + 0.8f * t), 0.f, 1.f)
                    };
                    const Vec3 a = trail[k - 1];
                    const Vec3 b = trail[k];
                    lines.push_back({a.x, a.y, a.z, b.x, b.y, b.z, tc.x, tc.y, tc.z});
                }
            }
        }
        for (auto it = trails.begin(); it != trails.end();)
        {
            const bool stillAlive = std::any_of(sim.organisms.begin(), sim.organisms.end(), [&](const Organism& o) { return o.id == it->first; });
            if (!stillAlive) it = trails.erase(it);
            else ++it;
        }

        // Environment layer: terrain + ocean + field probes
        if (showEnv)
        {
            for (int ix = -36; ix <= 36; ++ix)
            {
                for (int iz = -36; iz <= 36; ++iz)
                {
                    const float x = static_cast<float>(ix) * 3.8f;
                    const float z = static_cast<float>(iz) * 3.8f;
                    const Vec3 query{x, 0.f, z};
                    const float groundY = sim.env.groundHeight(query);
                    const float waterY = sim.env.waterSurfaceHeight(query);
                    const float shore = sim.env.shorelineBlend(Vec3{x, waterY, z});

                    const float slope = std::abs(sim.env.groundHeight(Vec3{x + 0.6f, 0.f, z})
                        - sim.env.groundHeight(Vec3{x - 0.6f, 0.f, z}))
                        + std::abs(sim.env.groundHeight(Vec3{x, 0.f, z + 0.6f})
                            - sim.env.groundHeight(Vec3{x, 0.f, z - 0.6f}));

                    const Biome biome = sim.env.biomeAt(Vec3{x, groundY + 0.6f, z});
                    const Vec3 biomeCol = biomeColor(biome);
                    const Vec3 landColor{
                        clampf(biomeCol.x + 0.18f * slope, 0.f, 1.f),
                        clampf(biomeCol.y + 0.20f * (1.f - slope * 0.2f) + 0.10f * (1.f - shore), 0.f, 1.f),
                        clampf(biomeCol.z + 0.08f * (1.f - shore), 0.f, 1.f)
                    };
                    pts.push_back({x, groundY, z, landColor.x, landColor.y, landColor.z, 2.1f});

                    if (waterY > groundY - 1.2f)
                    {
                        const float foam = clampf((waterY - groundY) / 4.0f, 0.f, 1.f);
                        const Vec3 waterColor{
                            clampf(0.06f + 0.18f * foam, 0.f, 1.f),
                            clampf(0.28f + 0.30f * shore, 0.f, 1.f),
                            clampf(0.45f + 0.42f * shore, 0.f, 1.f)
                        };
                        pts.push_back({x, waterY, z, waterColor.x, waterColor.y, waterColor.z, 2.4f});
                    }
                }
            }

            for (int i = 0; i < 340; i++)
            {
                const Vec3 h = hashToVec3(static_cast<uint64_t>(4100 + i));
                const Vec3 p = h * sim.worldRadius;
                if (len(p) > sim.worldRadius) continue;

                const float n = sim.env.nutrient(p);
                const float l = sim.env.light(p);
                const float x = sim.env.toxin(p);
                const float b = sim.env.biomass(p);

                const Vec3 col{
                    clampf(0.15f + 0.55f * x, 0.f, 1.f),
                    clampf(0.12f + 0.55f * n + 0.25f * b, 0.f, 1.f),
                    clampf(0.12f + 0.60f * l, 0.f, 1.f)
                };

                pts.push_back({p.x, p.y, p.z, col.x, col.y, col.z, 1.6f});
            }

            if (showCurrents)
            {
                for (int ix = -14; ix <= 14; ix += 2)
                {
                    for (int iy = -3; iy <= 4; iy += 2)
                    {
                        for (int iz = -14; iz <= 14; iz += 2)
                        {
                            Vec3 p{static_cast<float>(ix) * 6.f, static_cast<float>(iy) * 6.f, static_cast<float>(iz) * 6.f};
                            if (len(p) > sim.worldRadius * 0.92f) continue;

                            const Vec3 velocity = sim.env.fluidVelocity(p) + sim.env.windVelocity(p) * 0.45f;
                            const float speed = len(velocity);
                            if (speed < 0.03f) continue;

                            const Vec3 dir = velocity / speed;
                            const float arrowLen = clampf(1.5f + speed * 2.2f, 1.3f, 5.2f);
                            const Vec3 tip = p + dir * arrowLen;

                            const Vec3 col{
                                clampf(0.12f + 0.22f * speed, 0.f, 1.f),
                                clampf(0.55f + 0.20f * speed, 0.f, 1.f),
                                clampf(0.78f + 0.14f * speed, 0.f, 1.f)
                            };

                            lines.push_back({p.x, p.y, p.z, tip.x, tip.y, tip.z, col.x, col.y, col.z});
                        }
                    }
                }
            }

            // deterministic volumetric particles to make the ecosystem feel dense
            for (int ix = -8; ix <= 8; ++ix)
                for (int iy = -4; iy <= 5; ++iy)
                    for (int iz = -8; iz <= 8; ++iz)
                    {
                        const Vec3 seed = Vec3{static_cast<float>(ix), static_cast<float>(iy), static_cast<float>(iz)};
                        const Vec3 jitter = hashToVec3(static_cast<uint64_t>((ix + 19) * 73856093 ^ (iy + 23) * 19349663 ^ (iz + 31) * 83492791));
                        Vec3 p = seed * 10.0f + jitter * 3.0f;
                        p.y += std::sin(sim.env.time * 0.6f + seed.x * 0.37f + seed.z * 0.21f) * 2.5f;
                        if (len(p) > sim.worldRadius * 1.05f) continue;
                        const float n = sim.env.nutrient(p);
                        const float x = sim.env.toxin(p);
                        const float lum = clampf(0.04f + 0.13f * n + 0.1f * x, 0.f, 0.35f);
                        pts.push_back({p.x, p.y, p.z, lum * 0.8f, lum, lum * 1.1f, 1.2f + 1.4f * lum});
                    }
        }

        // celestial cues (sun/moon + stars) for time and atmosphere
        const float dayPhase = 0.5f + 0.5f * std::sin(sim.env.time * 0.07f);
        const float orbit = sim.env.time * 0.11f;
        const Vec3 sunPos{std::cos(orbit) * sim.worldRadius * 0.95f, 25.f + std::sin(orbit) * sim.worldRadius * 0.42f, std::sin(orbit) * sim.worldRadius * 0.95f};
        const Vec3 moonPos = sunPos * -0.95f + Vec3{0.f, 18.f, 0.f};
        pts.push_back({sunPos.x, sunPos.y, sunPos.z, 1.f, 0.84f, 0.4f, 16.f + 6.f * dayPhase});
        pts.push_back({moonPos.x, moonPos.y, moonPos.z, 0.55f, 0.68f, 1.f, 10.f + 5.f * (1.f - dayPhase)});
        for (int i = 0; i < 140; ++i)
        {
            const Vec3 h = hashToVec3(static_cast<uint64_t>(7000 + i));
            const Vec3 star = normalize(h) * (sim.worldRadius * 1.2f);
            const float twinkle = 0.65f + 0.35f * std::sin(sim.env.time * 1.3f + i * 3.1f);
            const float vis = (1.f - dayPhase) * twinkle;
            pts.push_back({star.x, star.y, star.z, 0.45f * vis, 0.55f * vis, 0.75f * vis, 1.5f + 1.8f * vis});
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
        const float dayL = 0.18f + 0.82f * dayPhase;
        glClearColor(0.01f + dayL * 0.06f, 0.02f + dayL * 0.1f, 0.05f + dayL * 0.12f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const Vec3 eye{
            camDist * std::cos(camPitch) * std::cos(camYaw),
            camDist * std::sin(camPitch),
            camDist * std::cos(camPitch) * std::sin(camYaw)
        };

        const Mat4 P = perspective(60.f * 3.1415926f / 180.f,
                                   static_cast<float>(w) / static_cast<float>(std::max(h, 1)),
                                   0.1f, 1000.f);
        Vec3 focus{0, 0, 0};
        if (cinematicCam && !sim.organisms.empty())
        {
            const auto target = std::max_element(sim.organisms.begin(), sim.organisms.end(), [](const Organism& a, const Organism& b) { return a.energy < b.energy; });
            focus = Sim::organismCenterOfMass(*target);
        }
        const Mat4 V  = lookAt(eye + focus, focus, Vec3{0, 1, 0});
        const Mat4 VP = mul(P, V);

        // Lines (springs)
        glUseProgram(progLines);
        glUniformMatrix4fv(uVPLines, 1, GL_FALSE, VP.m);
        glUniform1f(uFogLines, fogDensity);
        glBindVertexArray(vaoLines);
        glLineWidth(1.0f);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(lineVerts.size()));

        // Points (cells + env probes)
        glUseProgram(progPoints);
        glUniformMatrix4fv(uVPPoints, 1, GL_FALSE, VP.m);
        glUniform1f(uTimePoints, static_cast<float>(sim.env.time));
        glUniform1f(uFogPoints, fogDensity);
        glUniform1f(uGlowPoints, 0.8f + 0.4f * (1.f - dayPhase));
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

            std::ostringstream title;
            title << "EvoLife3D Cinematic | Pop " << sim.organisms.size()
                << " | Fitness " << std::fixed << std::setprecision(2) << sim.metrics.meanFitness
                << " | Het " << sim.metrics.heterozygosity
                << " | Div " << std::setprecision(2) << sim.metrics.nicheDiversity
                << " | Cells " << std::setprecision(1) << sim.metrics.meanCellCount
                << " | x" << std::setprecision(1) << timeScale
                << (cinematicCam ? " | CAM:follow" : " | CAM:free")
                << (showTrails ? " | trails" : " | no-trails");
            glfwSetWindowTitle(window, title.str().c_str());
        }

        if (appScreen == AppScreen::MainMenu)
        {
            std::ostringstream title;
            title << "EvoLife3D | MENU  [Enter] lancer  [O] options  [Esc] quitter"
                  << " | Population initiale: " << settings.initialPopulation
                  << " | Vitesse: x" << std::fixed << std::setprecision(1) << settings.startTimeScale;
            glfwSetWindowTitle(window, title.str().c_str());
        }
        else if (appScreen == AppScreen::Options)
        {
            std::ostringstream title;
            title << "EvoLife3D | OPTIONS  [I/K] population +/-  [1/2] vitesse +/-  [P] demarrer pause "
                  << (settings.startPaused ? "on" : "off")
                  << "  [G] env " << (settings.showEnv ? "on" : "off")
                  << "  [C] courants " << (settings.showCurrents ? "on" : "off")
                  << "  [F] camera " << (settings.cinematicCam ? "cinema" : "libre")
                  << "  [T] traces " << (settings.showTrails ? "on" : "off")
                  << "  [O] retour"
                  << " | Pop=" << settings.initialPopulation
                  << " | x" << std::fixed << std::setprecision(1) << settings.startTimeScale;
            glfwSetWindowTitle(window, title.str().c_str());
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
