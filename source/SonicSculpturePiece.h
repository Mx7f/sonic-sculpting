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

    CFrame          m_frame;

    void getTransformedFramesFromOriginals();
    void uploadToGPU();
public:
    
    SonicSculpturePiece() : m_gpuUpdated(false) {}

    void insert(const CFrame& frame, const float radius, const float delta);

    void setShaderArgs(Args& args);
    


};


#endif