#include "evo/Renderer.hpp"

#include <glad/glad.h>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace evo
{
    Renderer::~Renderer()
    {
        shutdown();
    }

    uint32_t Renderer::compileShader(const uint32_t type, const char *src)
    {
        const GLuint s = glCreateShader(type);
        glShaderSource(s, 1, &src, nullptr);
        glCompileShader(s);

        GLint ok = 0;
        glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
        if (!ok)
        {
            char log[4096];
            glGetShaderInfoLog(s, static_cast<GLsizei>(sizeof(log)), nullptr, log);
            std::cerr << "Shader compile error:\n" << log << "\n";
        }
        return s;
    }

    uint32_t Renderer::makeProgram(const char *vs, const char *fs)
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
            char log[4096];
            glGetProgramInfoLog(p, static_cast<GLsizei>(sizeof(log)), nullptr, log);
            std::cerr << "Program link error:\n" << log << "\n";
        }

        glDeleteShader(v);
        glDeleteShader(f);
        return p;
    }

    void Renderer::buildSphereMesh(int stacks, int slices)
    {
        // Simple UV-sphere. Low poly is fine for instancing.
        struct V
        {
            float px, py, pz;
            float nx, ny, nz;
        };

        std::vector<V>        verts;
        std::vector<uint32_t> idx;

        stacks = std::max(3, stacks);
        slices = std::max(3, slices);

        for (int i = 0; i <= stacks; ++i)
        {
            const float v   = static_cast<float>(i) / static_cast<float>(stacks);
            const float phi = v * 3.1415926535f; // 0..pi

            for (int j = 0; j <= slices; ++j)
            {
                const float u     = static_cast<float>(j) / static_cast<float>(slices);
                const float theta = u * 2.0f * 3.1415926535f; // 0..2pi

                const float x = std::sin(phi) * std::cos(theta);
                const float y = std::cos(phi);
                const float z = std::sin(phi) * std::sin(theta);

                verts.push_back({x, y, z, x, y, z});
            }
        }

        const int stride = slices + 1;
        for (int i = 0; i < stacks; ++i)
        {
            for (int j = 0; j < slices; ++j)
            {
                const auto a = static_cast<uint32_t>(i * stride + j);
                const auto b = static_cast<uint32_t>((i + 1) * stride + j);
                const auto c = static_cast<uint32_t>((i + 1) * stride + (j + 1));
                const auto d = static_cast<uint32_t>(i * stride + (j + 1));

                idx.push_back(a);
                idx.push_back(b);
                idx.push_back(c);

                idx.push_back(a);
                idx.push_back(c);
                idx.push_back(d);
            }
        }

        m_indexCount = static_cast<int>(idx.size());

        // GL buffers
        glGenVertexArrays(1, &m_vaoSphere);
        glGenBuffers(1, &m_vboSphere);
        glGenBuffers(1, &m_eboSphere);

        glBindVertexArray(m_vaoSphere);

        glBindBuffer(GL_ARRAY_BUFFER, m_vboSphere);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(verts.size() * sizeof(V)), verts.data(),
                     GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_eboSphere);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(idx.size() * sizeof(uint32_t)),
                     idx.data(),
                     GL_STATIC_DRAW);

        // attrib 0: pos
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(V), static_cast<void *>(nullptr));

        // attrib 1: normal
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(V),
                              reinterpret_cast<void *>(3 * sizeof(float)));

        // instance VBO will be bound later
        glBindVertexArray(0);
    }

    bool Renderer::init()
    {
        if (m_inited)
            return true;

        // Sphere instanced shader
        const auto VS_SPHERE = R"(
        #version 330 core
        layout(location=0) in vec3 aPos;
        layout(location=1) in vec3 aNor;

        // Instances:
        layout(location=2) in vec3 iPos;
        layout(location=3) in float iRadius;
        layout(location=4) in vec3 iColor;

        uniform mat4 uVP;

        out vec3 vColor;
        out vec3 vNor;

        void main(){
            vec3 wp = iPos + aPos * iRadius;
            vNor = aNor; // sphere normals in object space == world space here (uniform scale)
            vColor = iColor;
            gl_Position = uVP * vec4(wp, 1.0);
        }
    )";

        const auto FS_SPHERE = R"(
        #version 330 core
        in vec3 vColor;
        in vec3 vNor;

        uniform vec3 uLightDir; // should be normalized

        out vec4 FragColor;

        void main(){
            float ndl = max(dot(normalize(vNor), normalize(uLightDir)), 0.0);
            float amb = 0.22;
            float lit = amb + (1.0-amb) * ndl;

            FragColor = vec4(vColor * lit, 1.0);
        }
    )";

        // Line shader (debug)
        const auto VS_LINE = R"(
        #version 330 core
        layout(location=0) in vec3 aPos;
        layout(location=1) in vec3 aColor;

        uniform mat4 uVP;
        out vec3 vColor;

        void main(){
            gl_Position = uVP * vec4(aPos, 1.0);
            vColor = aColor;
        }
    )";

        const auto FS_LINE = R"(
        #version 330 core
        in vec3 vColor;
        out vec4 FragColor;

        void main(){
            FragColor = vec4(vColor, 1.0);
        }
    )";

        m_progSphere = makeProgram(VS_SPHERE, FS_SPHERE);
        m_uVP_sphere = glGetUniformLocation(m_progSphere, "uVP");
        m_uLightDir  = glGetUniformLocation(m_progSphere, "uLightDir");

        m_progLine = makeProgram(VS_LINE, FS_LINE);
        m_uVP_line = glGetUniformLocation(m_progLine, "uVP");

        buildSphereMesh(10, 10);

        // Instance buffer
        glGenBuffers(1, &m_vboInstances);
        glBindVertexArray(m_vaoSphere);
        glBindBuffer(GL_ARRAY_BUFFER, m_vboInstances);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sizeof(RenderPoint) * 1024), nullptr,
                     GL_DYNAMIC_DRAW);

        // Instance attribs:
        // location 2: iPos (vec3)
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(RenderPoint), static_cast<void *>(nullptr));
        glVertexAttribDivisor(2, 1);

        // location 3: iRadius (float) offset 6 floats
        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(RenderPoint),
                              reinterpret_cast<void *>(6 * sizeof(float)));
        glVertexAttribDivisor(3, 1);

        // location 4: iColor (vec3) offset 3 floats
        glEnableVertexAttribArray(4);
        glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(RenderPoint),
                              reinterpret_cast<void *>(3 * sizeof(float)));
        glVertexAttribDivisor(4, 1);

        glBindVertexArray(0);

        // Lines VAO/VBO
        glGenVertexArrays(1, &m_vaoLines);
        glGenBuffers(1, &m_vboLines);
        glBindVertexArray(m_vaoLines);
        glBindBuffer(GL_ARRAY_BUFFER, m_vboLines);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sizeof(LineVertex) * 1024), nullptr,
                     GL_DYNAMIC_DRAW);

        // pos
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(LineVertex), static_cast<void *>(nullptr));
        // color
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(LineVertex),
                              reinterpret_cast<void *>(3 * sizeof(float)));

        glBindVertexArray(0);

        m_inited = true;
        return true;
    }

    void Renderer::shutdown()
    {
        if (!m_inited)
            return;

        if (m_vboLines)
            glDeleteBuffers(1, &m_vboLines);
        if (m_vaoLines)
            glDeleteVertexArrays(1, &m_vaoLines);

        if (m_vboInstances)
            glDeleteBuffers(1, &m_vboInstances);

        if (m_eboSphere)
            glDeleteBuffers(1, &m_eboSphere);
        if (m_vboSphere)
            glDeleteBuffers(1, &m_vboSphere);
        if (m_vaoSphere)
            glDeleteVertexArrays(1, &m_vaoSphere);

        if (m_progLine)
            glDeleteProgram(m_progLine);
        if (m_progSphere)
            glDeleteProgram(m_progSphere);

        m_progLine   = 0;
        m_progSphere = 0;

        m_inited = false;
    }

    void Renderer::resize(const int w, const int h)
    {
        m_w = std::max(1, w);
        m_h = std::max(1, h);
        glViewport(0, 0, m_w, m_h);
    }

    void Renderer::setCameraVP(const Mat4 &vp)
    {
        m_vp = vp;
    }

    void Renderer::uploadInstances(const std::vector<RenderPoint> &instances)
    {
        m_instanceCount = static_cast<int>(instances.size());
        glBindBuffer(GL_ARRAY_BUFFER, m_vboInstances);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(instances.size() * sizeof(RenderPoint)),
                     instances.data(),
                     GL_DYNAMIC_DRAW);
    }

    void Renderer::drawSpheres() const
    {
        if (m_instanceCount <= 0)
            return;

        glUseProgram(m_progSphere);
        glUniformMatrix4fv(m_uVP_sphere, 1, GL_FALSE, m_vp.m);

        // A fixed light direction (camera-ish). Keep it simple.
        glUniform3f(m_uLightDir, -0.35f, 0.8f, 0.35f);

        glBindVertexArray(m_vaoSphere);
        glDrawElementsInstanced(GL_TRIANGLES, static_cast<GLsizei>(m_indexCount), GL_UNSIGNED_INT, nullptr,
                                static_cast<GLsizei>(m_instanceCount));
        glBindVertexArray(0);
    }

    void Renderer::uploadLines(const std::vector<LineVertex> &lines)
    {
        m_lineVertexCount = static_cast<int>(lines.size());
        glBindBuffer(GL_ARRAY_BUFFER, m_vboLines);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(lines.size() * sizeof(LineVertex)),
                     lines.data(),
                     GL_DYNAMIC_DRAW);
    }

    void Renderer::drawLines() const
    {
        if (!m_showLines)
            return;
        if (m_lineVertexCount <= 1)
            return;

        glUseProgram(m_progLine);
        glUniformMatrix4fv(m_uVP_line, 1, GL_FALSE, m_vp.m);

        glBindVertexArray(m_vaoLines);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(m_lineVertexCount));
        glBindVertexArray(0);
    }
} // namespace evo