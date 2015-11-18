#ifndef SonicSculpturePiece_h
#define SonicSculpturePiece_h
#include <G3D/G3DAll.h>

class SonicSculpturePiece {
protected:
    Array<float>    m_deltas;
    Array<CFrame>   m_originalFrames;
    Array<CFrame>   m_transformedFrames;
    Array<float>    m_radii;

    bool            m_gpuUpdated;

    AttributeArray  m_vertices;
    AttributeArray  m_texCoord;
    AttributeArray  m_normals;
    IndexStream     m_indices;

    shared_ptr<UniversalMaterial> m_material;

    CFrame          m_frame;

    void getTransformedFramesFromOriginals();
    void uploadToGPU();
    void getCPUGeometry(CPUVertexArray& cpuVertexArray, Array<int>& cpuIndexArray) const;
    SonicSculpturePiece() : m_gpuUpdated(false) {}
public:

    static shared_ptr<SonicSculpturePiece> create(shared_ptr<UniversalMaterial> material);

    void insert(const CFrame& frame, const float radius, const float delta);

    void render(RenderDevice* rd, const LightingEnvironment& env);

    void setShaderArgs(Args& args);

    int size() const {
        return m_originalFrames.size();
    }
    
};


#endif