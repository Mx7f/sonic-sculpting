#ifndef SonicSculpturePiece_h
#define SonicSculpturePiece_h
#include <G3D/G3DAll.h>
#include "AudioSample.h"
class SonicSculpturePiece {
protected:
    CFrame          m_frame;
    Array<float>	m_audioSamples;
    Array<float>    m_deltas;
    Array<CFrame>   m_originalFrames;

    Array<CFrame>   m_transformedFrames;
    Array<CFrame>   m_previousTransformedFrames;

    Array<float>    m_radii;

    Vector3         m_generalDirection;

    bool            m_gpuUpdated;

    AttributeArray  m_vertices;
    AttributeArray  m_texCoord;
    AttributeArray  m_normals;
    AttributeArray  m_boneIndices;
    AttributeArray  m_boneWeights;
    IndexStream     m_indices;

    shared_ptr<Texture> m_boneTexture;

    shared_ptr<UniversalMaterial> m_material;

    bool m_frozen;

    void uploadToGPU();
    void getCPUGeometry(CPUVertexArray& cpuVertexArray, Array<int>& cpuIndexArray) const;
    SonicSculpturePiece() : m_gpuUpdated(false) {}
public:
    shared_ptr<UniversalMaterial> material() const {
        return m_material;
    }

    void onSimulation(RealTime rdt, SimTime sdt, SimTime idt);

    void serialize(BinaryOutput& output) const;

    shared_ptr<AudioSample> getBaseAudioSample() const;

	shared_ptr<AudioSample> getAudioSampleFromRay(const Ray& ray);

    static shared_ptr<SonicSculpturePiece> create(shared_ptr<UniversalMaterial> material);

    void insert(const CFrame& frame, const float radius, const float delta, const Array<float>& samples = Array<float>());

    void render(RenderDevice* rd, const LightingEnvironment& env);

    void setShaderArgs(Args& args);

    static shared_ptr<SonicSculpturePiece> fromBinaryInput(shared_ptr<UniversalMaterial> material, BinaryInput& input);

    int size() const {
        return m_originalFrames.size();
    }

    void setFrozen(bool b) {
        m_frozen = b;
    }
    
};


#endif