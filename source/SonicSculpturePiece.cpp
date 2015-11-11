#include "SonicSculpturePiece.h"

void SonicSculpturePiece::getTransformedFramesFromOriginals() {
    m_transformedFrames.fastClear();
    m_transformedFrames.resize(m_originalFrames.size());
    float totalDelta = 0.0f;
    // TODO: Find more elegant algorithm
    for (int i = m_originalFrames.size() - 1; i >= 0; --i) {
        m_transformedFrames[i] = m_originalFrames[i];
        m_transformedFrames[i].translation += m_originalFrames[i].lookVector() * totalDelta;
        totalDelta += m_deltas[i];
    }
}

void SonicSculpturePiece::getCPUGeometry(CPUVertexArray & cpuVertexArray, Array<int>& cpuIndexArray) const {
    int sideCount = 6;
    int segmentCount = m_transformedFrames.size();
    int currentPointIndex = 0;
    cpuVertexArray.vertex.resize(segmentCount*sideCount);
    for (int j = 0; j < segmentCount; ++j) {
        const CFrame& frame = m_transformedFrames[j];
        const float thickness = m_radii[j];
        const Point3& currentPoint = frame.translation;

        for (int k = 0; k < sideCount; ++k) {
            CPUVertexArray::Vertex& v = cpuVertexArray.vertex[currentPointIndex * sideCount + k];

            // Rotate around the axis of direction
            v.normal = Matrix3::fromAxisAngle(-frame.lookVector(), k * 2.0f * pif() / sideCount) * frame.upVector();
            v.position = currentPoint + v.normal * thickness;
            v.texCoord0 = Point2(float(k) / sideCount, float(j) / float(segmentCount));
            v.tangent = Vector4(frame.lookVector(), 1.0f);
        }

        if (j > 0) {
            int currStart = currentPointIndex * sideCount;
            int prevStart = currStart - sideCount;
            for (int k = 0; k < sideCount; ++k) {
                const int next = (k + 1) % sideCount;
                cpuIndexArray.append(prevStart + k, prevStart + next, currStart + k);
                cpuIndexArray.append(prevStart + next, currStart + next, currStart + k);
            }
        }

        ++currentPointIndex;
    }
}

void SonicSculpturePiece::uploadToGPU() {
    getTransformedFramesFromOriginals();

    CPUVertexArray cpuVertexArray;
    Array<int> cpuIndexArray;
    getCPUGeometry(cpuVertexArray, cpuIndexArray);
    cpuVertexArray.copyToGPU(m_vertices, m_normals, AttributeArray(), m_texCoord, AttributeArray(), AttributeArray(), VertexBuffer::WRITE_ONCE);
    m_indices = IndexStream(cpuIndexArray, VertexBuffer::create(sizeof(int)*cpuIndexArray.size() + 8));
    
    m_gpuUpdated = true;
}

shared_ptr<SonicSculpturePiece> SonicSculpturePiece::create(shared_ptr<UniversalMaterial> material) {
    shared_ptr<SonicSculpturePiece> s(new SonicSculpturePiece());
    s->m_material = material;
    return s;
}

void SonicSculpturePiece::insert(const CFrame & frame, const float radius, const float delta) {
    m_gpuUpdated = false;
    m_originalFrames.append(frame);
    m_radii.append(radius);
    m_deltas.append(delta);
    uploadToGPU();
}

void SonicSculpturePiece::render(RenderDevice * rd, const LightingEnvironment& env) {
    if (m_radii.size() < 2) {
        return;
    }
    Args args;
    env.setShaderArgs(args);
    rd->setObjectToWorldMatrix(m_frame);
    args.setMacro("OPAQUE_PASS", 1);
    setShaderArgs(args);
    rd->setCullFace(CullFace::NONE);
    LAUNCH_SHADER("UniversalSurface/UniversalSurface_render.*", args);
}

void SonicSculpturePiece::setShaderArgs(Args & args) {
    if (!m_gpuUpdated) {
        uploadToGPU();
    }
    args.setAttributeArray("g3d_Vertex", m_vertices);
    args.setAttributeArray("g3d_Normal", m_normals);
    args.setAttributeArray("g3d_TexCoord0", m_texCoord);
    args.setIndexStream(m_indices);
    args.setMacro("NUM_BONES", 0);
    args.setMacro("NUM_LIGHTMAP_DIRECTIONS", 0);
    args.setMacro("HAS_VERTEX_COLOR", 0);
    args.setMacro("HAS_ALPHA", m_material->hasAlpha());
    args.setMacro("ALPHA_HINT", m_material->alphaHint());
    args.setMacro("HAS_TRANSMISSIVE", m_material->hasTransmissive());
    args.setMacro("HAS_EMISSIVE", false);

    m_material->setShaderArgs(args, "material.");
}
