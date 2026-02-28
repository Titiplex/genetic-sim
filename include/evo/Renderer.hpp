#pragma once
#include <cstdint>
#include <vector>

#include "Math.hpp"
#include "RenderTypes.hpp"

struct GLFWwindow;

namespace evo
{
    class Renderer
    {
        public:
            Renderer() = default;
            ~Renderer();

            bool init();
            void shutdown();

            void resize(int w, int h);
            void setCameraVP(const Mat4 &vp);

            // Instanced spheres: one instance per cell
            void uploadInstances(const std::vector<RenderPoint> &instances);
            void drawSpheres() const;

            // Debug: springs/lines
            void uploadLines(const std::vector<LineVertex> &lines); // pairs (A,B)
            void drawLines() const;

            // Simple toggles
            void setShowLines(const bool v)
            {
                m_showLines = v;
            }

            [[nodiscard]] bool showLines() const
            {
                return m_showLines;
            }

        private:
            bool m_inited = false;

            int  m_w  = 1;
            int  m_h  = 1;
            Mat4 m_vp = Mat4::identity();

            bool m_showLines = false;

            // GL programs
            uint32_t m_progSphere = 0;
            int      m_uVP_sphere = -1;
            int      m_uLightDir  = -1;

            uint32_t m_progLine = 0;
            int      m_uVP_line = -1;

            // Sphere mesh
            uint32_t m_vaoSphere  = 0;
            uint32_t m_vboSphere  = 0;
            uint32_t m_eboSphere  = 0;
            int      m_indexCount = 0;

            // Instance buffer
            uint32_t m_vboInstances  = 0;
            int      m_instanceCount = 0;

            // Lines
            uint32_t m_vaoLines        = 0;
            uint32_t m_vboLines        = 0;
            int      m_lineVertexCount = 0;

        private:
            static uint32_t compileShader(uint32_t type, const char *src);
            static uint32_t makeProgram(const char *vs, const char *fs);

            void buildSphereMesh(int stacks, int slices);
    };
} // namespace evo